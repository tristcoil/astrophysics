"""The Moon crossing Earth's face: three transits, and how bright the Moon is.

EPIC's gallery holds the natural-colour frames of the lunar transits of
16 July 2015, 5 July 2016 and 11 February 2021 (the last also in the regular
archive, pixel for pixel the same). DSCOVR sits sunward of Earth, so in a
transit it sees the Moon's fully lit far side, the hemisphere never seen from
Earth.

For every frame we predict where the Moon must appear from JPL Horizons
vectors (DSCOVR, Moon, Sun) and the images' fixed plate scale, then measure
the Moon's inner disc (inside 0.8 of its radius: the three colour bands are
exposed minutes apart while the Moon moves ~40 pixels, which fringes its
edges) and Earth's inner disc away from the Moon.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

import ephemeris as E
import geometry as G
import photometry as ph
from config import GALLERY_DIR, PLATE_SCALE_1080
from transform import GAMMA_GRID, GAMMA_NOMINAL, YANG2018

TRANSITS = {
    "2015": "lunar_transit_2015",
    "2016": "lunar_transit_2016",
    "2021": "lunar_transit_2021",
}
MOON_INNER = 0.8  # fraction of the Moon's radius used
EARTH_INNER = 0.978  # as for the Earth photometry (solar zenith < 78 deg)
MOON_MARGIN_PX = 60  # Earth pixels this close to the Moon's limb are skipped (2048-px frame)

#: EPIC's own whole-disc average reflectance of the lunar far side, read off
#: Gorkavyi, Krotkov & Marshak 2026 (Front. Remote Sens. 7, 1756531), Fig. 7,
#: curve 3, from images of the Moon against space (+-0.003 from the reading).
FARSIDE_REFLECTANCE = {"blue": 0.105, "green": 0.130, "red": 0.151}


def frame_time(name: str) -> pd.Timestamp:
    """Time from a gallery file name: epic_1b_YYYYMMDDhhmmss_vv.png or the 2015
    style DOY_YYYYDOYhhmmss.png."""
    m = re.search(r"epic_1b_(\d{14})", name)
    if m:
        return pd.to_datetime(m.group(1), format="%Y%m%d%H%M%S")
    m = re.search(r"\d{3}_(\d{4})(\d{3})(\d{6})", name)
    year, doy, hms = m.groups()
    return (pd.Timestamp(f"{year}-01-01") + pd.Timedelta(days=int(doy) - 1)
            + pd.to_timedelta(f"{hms[:2]}:{hms[2:4]}:{hms[4:]}"))


def frames(key: str) -> pd.DataFrame:
    files = sorted((GALLERY_DIR / TRANSITS[key]).glob("*.png"))
    return pd.DataFrame({"path": files, "t": [frame_time(p.name) for p in files]}).sort_values("t")


def geometry_at(times: pd.Series) -> pd.DataFrame:
    """Horizons geometry interpolated to the frame times (2-minute grid)."""
    t0 = (times.min() - pd.Timedelta("10min")).strftime("%Y-%m-%d %H:%M")
    t1 = (times.max() + pd.Timedelta("10min")).strftime("%Y-%m-%d %H:%M")
    vec = E.vectors(t0, t1, E.FINE_STEP)
    x = (vec["t"] - vec["t"].iloc[0]).dt.total_seconds().to_numpy()
    xi = (times - vec["t"].iloc[0]).dt.total_seconds().to_numpy()
    out = pd.DataFrame({"t": times.to_numpy()})
    for col in vec.columns.drop("t"):
        out[col] = np.interp(xi, x, vec[col].to_numpy())
    return out


def fit_earth_limb(rgb: np.ndarray, exclude_xy=None, exclude_r: float = 0.0):
    """Centre and radius of Earth's disc from its limb (algebraic circle fit).

    Limb points are disc pixels next to masked space, away from the frame
    edges and from the Moon (which may stick out past Earth's limb)."""
    mask = ph.disc_mask(rgb)
    n = mask.shape[0]
    edge = mask & ~(np.roll(mask, 1, 0) & np.roll(mask, -1, 0) & np.roll(mask, 1, 1) & np.roll(mask, -1, 1))
    yy, xx = np.nonzero(edge)
    keep = (xx > 5) & (xx < n - 6) & (yy > 5) & (yy < n - 6)
    if exclude_xy is not None:
        keep &= np.hypot(xx - exclude_xy[0], yy - exclude_xy[1]) > exclude_r
    x, y = xx[keep].astype(float), yy[keep].astype(float)
    a = np.column_stack([x, y, np.ones_like(x)])
    b = x**2 + y**2
    sol = np.linalg.lstsq(a, b, rcond=None)[0]
    cx, cy = sol[0] / 2, sol[1] / 2
    r = np.sqrt(sol[2] + cx**2 + cy**2)
    # one robust pass: drop points far from the first circle
    d = np.abs(np.hypot(x - cx, y - cy) - r)
    k = d < 3 * np.median(d) + 1
    sol = np.linalg.lstsq(a[k], b[k], rcond=None)[0]
    cx, cy = sol[0] / 2, sol[1] / 2
    return float(cx), float(cy), float(np.sqrt(sol[2] + cx**2 + cy**2))


def moon_likeness(rgb: np.ndarray) -> np.ndarray:
    """Pixels coloured like the Moon in these images: grey-brown, a little redder
    than blue (red - blue 6..32), moderately dark (red 35..110). Ocean is bluer,
    clouds greyer and brighter; desert can pass, which the disc filter handles."""
    r, g, b = (rgb[..., i].astype(np.int16) for i in range(3))
    return ((r - b >= 6) & (r - b <= 32) & (r >= 35) & (r <= 110) & (g >= 30)).astype(np.float32)


def locate_moon(rgb: np.ndarray, guess_xy, r_moon_px: float, search_px: float = 700.0):
    """Find the Moon's centre: correlate the Moon-likeness map with a disc of
    the Moon's size near the predicted position, then refine by centroiding."""
    from scipy.signal import fftconvolve

    like = moon_likeness(rgb)
    step = 4  # work on a 4x-reduced map for speed
    small = like[::step, ::step]
    rk = r_moon_px / step
    yk, xk = np.mgrid[-int(rk) - 1:int(rk) + 2, -int(rk) - 1:int(rk) + 2]
    kern = (np.hypot(xk, yk) <= 0.85 * rk).astype(np.float32)
    kern /= kern.sum()
    score = fftconvolve(small, kern, mode="same")
    ny, nx = score.shape
    gy, gx = np.mgrid[:ny, :nx]
    near = np.hypot(gx * step - guess_xy[0], gy * step - guess_xy[1]) < search_px
    score = np.where(near, score, -1)
    j = np.unravel_index(np.argmax(score), score.shape)
    cx, cy = j[1] * step, j[0] * step
    yy, xx = np.mgrid[:rgb.shape[0], :rgb.shape[1]]
    for _ in range(3):
        sel = (np.hypot(xx - cx, yy - cy) < 0.9 * r_moon_px) & (like > 0)
        if sel.sum() < 100:
            break
        cx, cy = xx[sel].mean(), yy[sel].mean()
    return float(cx), float(cy), float(score.max())


def measure_transit(key: str) -> pd.DataFrame:
    fr = frames(key)
    geo = geometry_at(fr["t"])
    dsc = geo[["dscovr_x", "dscovr_y", "dscovr_z"]].to_numpy()
    moon = geo[["moon_x", "moon_y", "moon_z"]].to_numpy()
    sun = geo[["sun_x", "sun_y", "sun_z"]].to_numpy()
    mv = G.moon_in_view(dsc, moon, sun)
    v = np.arange(256) / 255.0
    rows = []
    for k, (_, f) in enumerate(fr.iterrows()):
        rgb = ph.load_rgb(f["path"])
        n = rgb.shape[0]
        scale = PLATE_SCALE_1080 * n / 1080.0
        c = (n - 1) / 2.0
        re_pred = mv["r_earth_rad"][k] * scale
        rm_pred = mv["r_moon_rad"][k] * scale
        px, py = c + mv["x_rad"][k] * scale, c - mv["y_rad"][k] * scale
        # Locate the Moon in the frame (the 2015 set is not north-up and not
        # centred; the 2016 set is off by tens of pixels), then fit Earth's limb.
        mx, my, _ = locate_moon(rgb, (px, py), rm_pred, search_px=900.0)
        ex, ey, re_px = fit_earth_limb(rgb, (mx, my), 1.25 * rm_pred)
        rm_px = rm_pred * re_px / re_pred  # same plate scale for Moon and Earth
        yy, xx = np.mgrid[:n, :n]
        de = np.hypot(xx - ex, yy - ey)
        dm = np.hypot(xx - mx, yy - my)
        moon_px = dm < MOON_INNER * rm_px
        on_earth = moon_px & (de < re_px - 20)
        on_space = moon_px & (de > re_px + 20)
        earth_px = (de < EARTH_INNER * re_px) & (dm > rm_px + MOON_MARGIN_PX * n / 2048)
        # Region A: a fixed patch of lunar terrain (east part of the disc), to
        # compare the same ground in front of Earth and in front of space.
        region_a = moon_px & (xx - mx > 0.35 * rm_px) & (np.abs(yy - my) < 0.45 * rm_px)
        # Where the ephemeris puts the Moon relative to Earth's fitted centre,
        # vs where it is: tests the images' orientation and the geometry.
        sep_pred = np.hypot(mv["x_rad"][k], mv["y_rad"][k]) * scale * re_px / re_pred
        sep_meas = np.hypot(mx - ex, my - ey)
        row = {"transit": key, "file": f["path"].name, "t": f["t"],
               "moon_x_px": mx, "moon_y_px": my, "moon_r_px": rm_px, "earth_r_px": re_px,
               "earth_x_px": ex, "earth_y_px": ey, "earth_r_pred_px": re_pred,
               "moon_x_pred_px": px, "moon_y_pred_px": py,
               "sep_pred_px": sep_pred, "sep_meas_px": sep_meas,
               "pa_pred_deg": float(np.degrees(np.arctan2(mv["x_rad"][k], mv["y_rad"][k]))),
               "pa_meas_deg": float(np.degrees(np.arctan2(mx - ex, ey - my))),
               "moon_phase_deg": mv["moon_phase_deg"][k], "earth_phase_deg":
               float(G.phase_angle_deg(dsc[k:k + 1], sun[k:k + 1])[0]),
               "n_moon_on_earth": int(on_earth.sum()), "n_moon_on_space": int(on_space.sum()),
               "n_earth": int(earth_px.sum()),
               "regionA_on_space": float((region_a & (de > re_px + 20)).sum() / max(region_a.sum(), 1)),
               "regionA_on_earth": float((region_a & (de < re_px - 20)).sum() / max(region_a.sum(), 1))}
        for region, m in (("moon", moon_px), ("moonE", on_earth), ("moonS", on_space),
                          ("earth", earth_px), ("regionA", region_a)):
            if m.sum() < 2000:
                continue
            for ci, ch in enumerate(ph.CHANNELS):
                hist = np.bincount(rgb[..., ci][m], minlength=256)
                for g in GAMMA_GRID:
                    row[f"{region}_{ch}_g{g:.1f}"] = float(hist @ (v ** g) / hist.sum())
        rows.append(row)
    return pd.DataFrame(rows)


def summarise(tab: pd.DataFrame, earth_july2016: dict) -> dict:
    """Moon/Earth brightness per transit and channel, at every stored exponent.

    ``earth_july2016[gamma][channel]``: Earth's mean linear light over Yang et
    al.'s 12 days (inner disc), which fixes the reflectance scale."""
    out = {}
    for key, t in tab.groupby("transit"):
        full = t[t["n_moon_on_space"] == 0]  # Moon wholly in front of Earth
        use = full if len(full) else t
        res = {"n_frames": int(len(t)), "n_full_on_earth": int(len(full)),
               "moon_phase_deg": float(t["moon_phase_deg"].median()),
               "earth_phase_deg": float(t["earth_phase_deg"].median()),
               "moon_diameter_px": float(2 * t["moon_r_px"].median()),
               "earth_diameter_px": float(2 * t["earth_r_px"].median())}
        for ch in ph.CHANNELS:
            d_moon = float(use[f"moon_{ch}_g1.0"].median()) * 255
            d_earth = float(use[f"earth_{ch}_g1.0"].median()) * 255
            res[f"display_{ch}"] = {"moon": d_moon, "earth": d_earth, "ratio": d_moon / d_earth}
            ratios, refl = {}, {}
            for g in GAMMA_GRID:
                lm = float(use[f"moon_{ch}_g{g:.1f}"].median())
                le = float(use[f"earth_{ch}_g{g:.1f}"].median())
                ratios[f"{g:.1f}"] = lm / le
                if ch in YANG2018["mean"]:
                    k = YANG2018["mean"][ch] / earth_july2016[f"{g:.1f}"][ch]
                    refl[f"{g:.1f}"] = k * lm
            res[f"moon_over_earth_{ch}"] = ratios
            if refl:
                res[f"moon_reflectance_{ch}"] = refl
        out[key] = res
    return out


def stray_light(tab: pd.DataFrame, key: str = "2021") -> dict:
    """The same lunar terrain (region A) in front of Earth and in front of space.

    Light scattered inside the camera from the bright Earth adds to whatever
    lies in front of it; Gorkavyi et al. (2026) left out Moon images taken
    against Earth for this reason. Region A is wholly in front of Earth early
    in the 2021 transit and wholly beyond Earth's limb at its end."""
    t = tab[tab["transit"] == key]
    on_e = t[t["regionA_on_earth"] > 0.99]
    on_s = t[t["regionA_on_space"] > 0.99]
    if len(on_e) == 0 or len(on_s) == 0:
        return {}
    out = {"n_frames_on_earth": int(len(on_e)), "n_frames_on_space": int(len(on_s))}
    for ch in ph.CHANNELS:
        a = on_e[f"regionA_{ch}_g1.0"].median() * 255
        b = on_s[f"regionA_{ch}_g1.0"].median() * 255
        la = on_e[f"regionA_{ch}_g{GAMMA_NOMINAL:.1f}"].median()
        lb = on_s[f"regionA_{ch}_g{GAMMA_NOMINAL:.1f}"].median()
        out[ch] = {"display_on_earth": float(a), "display_on_space": float(b),
                   "excess_display": float(a / b - 1), "excess_linear": float(la / lb - 1)}
    return out


def expected_ratios() -> dict:
    """Moon/Earth reflectance ratio expected from calibrated EPIC data: EPIC's
    far side (Gorkavyi et al. 2026) over Earth's July-2016 mean (Yang et al.
    2018). Green has no published Earth value; it is interpolated linearly in
    wavelength between 443 and 680 nm and flagged as such."""
    earth = dict(YANG2018["mean"])
    earth["green"] = earth["blue"] + (earth["red"] - earth["blue"]) * (551 - 443) / (680 - 443)
    return {ch: {"moon": FARSIDE_REFLECTANCE[ch], "earth": earth[ch],
                 "ratio": FARSIDE_REFLECTANCE[ch] / earth[ch],
                 "earth_interpolated": ch == "green"} for ch in ph.CHANNELS}


def gamma_matching_moon(summary: dict, key: str, ch: str) -> float | None:
    """The exponent at which our Moon reflectance equals EPIC's far-side value."""
    refl = summary[key].get(f"moon_reflectance_{ch}")
    if not refl:
        return None
    g = np.array([float(k) for k in refl])
    r = np.array(list(refl.values()))
    target = FARSIDE_REFLECTANCE[ch]
    order = np.argsort(g)
    g, r = g[order], r[order]
    if not (r.min() <= target <= r.max()):
        return None
    return float(np.interp(target, r[::-1], g[::-1]))  # r decreases with gamma


__all__ = ["measure_transit", "summarise", "stray_light", "expected_ratios",
           "gamma_matching_moon", "TRANSITS", "GAMMA_NOMINAL"]
