"""Earth as a point of light: one row per image, with geometry and brightness.

For each image the "dot's" brightness in a colour channel is the disc's total
linear light divided by the area of the full projected disc,

    F = sum over disc pixels of (V/255)^gamma  /  (pi a b),

where a and b are Earth's projected equatorial and polar radii in pixels,
predicted from DSCOVR's distance and the fixed plate scale. F is what a
distant telescope would measure if it could not resolve the planet, in units
of the light of a uniformly bright disc of the same size (an "apparent
albedo" in display-linear units; see ``transform.py`` for the conversion to
reflectance). Dividing by the predicted area, not the lit pixel count, keeps
the dark crescent of the phase angle in the measurement, as it would be for a
real exoplanet.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import epicdata
import geometry as G
import photometry as ph
from config import (
    DISC_AREA_TOLERANCE, EXCLUDE_FROM_LIGHTCURVES, JUMP_NEIGHBOUR_H, MAX_CENTRE_OFFSET_PX,
    MAX_DISPLAY_JUMP, PLATE_SCALE_1080, R_EARTH_EQ_KM, R_EARTH_POL_KM,
)
from transform import GAMMA_GRID, GAMMA_NOMINAL

CHANNELS = ph.CHANNELS


def add_geometry(df: pd.DataFrame, size_px: int) -> pd.DataFrame:
    """Phase angle, distances, sub-solar point and the predicted disc."""
    dsc, moon, sun = epicdata.vectors(df)
    df = df.copy()
    d = G.norm(dsc)
    df["dist_km"] = d
    df["sun_au"] = G.sun_earth_distance_au(sun)
    df["phase_deg"] = G.phase_angle_deg(dsc, sun)
    gmst = G.gmst_rad(df["t"])
    df["subsolar_lat"], df["subsolar_lon"] = G.sub_point(sun, gmst)
    df["subobs_lat_calc"], df["subobs_lon_calc"] = G.sub_point(dsc, gmst)
    scale = PLATE_SCALE_1080 * size_px / 1080.0
    a = scale * np.arcsin(R_EARTH_EQ_KM / d)
    phi = np.radians(df["lat"].to_numpy())
    b_km = np.hypot(R_EARTH_POL_KM * np.cos(phi), R_EARTH_EQ_KM * np.sin(phi))
    b = scale * np.arcsin(b_km / d)
    df["disc_area_px"] = np.pi * a * b
    df["r_pred_px"] = np.sqrt(a * b)
    mv = G.moon_in_view(dsc, moon, sun)
    half = np.radians(0.31)
    df["moon_in_fov"] = ((np.abs(mv["x_rad"]) < half + mv["r_moon_rad"])
                         & (np.abs(mv["y_rad"]) < half + mv["r_moon_rad"]))
    return df


def add_fluxes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for g in GAMMA_GRID:
        for c in CHANNELS:
            total = df[ph.moment_column("disc", c, g)] * df["n_disc"]
            df[f"F_{c}_g{g:.1f}"] = total / df["disc_area_px"]
    for c in CHANNELS:
        df[f"F_{c}"] = df[f"F_{c}_g{GAMMA_NOMINAL:.1f}"]
        df[f"D_{c}"] = df[ph.moment_column("disc", c, 1.0)]  # display units
    # Colour: blue/red light ratio at the nominal exponent.
    df["blue_red"] = df["F_blue"] / df["F_red"]
    return df


def quality(df: pd.DataFrame) -> pd.DataFrame:
    """Flag images whose disc is off-centre or the wrong size, and special days."""
    df = df.copy()
    centre = (df["nx"] - 1) / 2.0
    k = 1080.0 / df["nx"]
    df["centre_offset_px"] = np.hypot(df["cx"] - centre, df["cy"] - centre) * k
    lit = df["n_disc"] / (df["disc_area_px"] * G.sunlit_fraction_of_disc(df["phase_deg"]))
    df["lit_area_ratio"] = lit
    special = df["t"].dt.strftime("%Y-%m-%d").isin(EXCLUDE_FROM_LIGHTCURVES)
    good = ((df["centre_offset_px"] < MAX_CENTRE_OFFSET_PX)
            & (np.abs(lit - 1.0) < DISC_AREA_TOLERANCE)
            & ~df["moon_in_fov"] & ~special)
    df["jump"] = frame_jump(df, good)
    df["ok"] = good & ~(np.abs(df["jump"]) > MAX_DISPLAY_JUMP)
    return df


def frame_jump(df: pd.DataFrame, good: pd.Series) -> pd.Series:
    """Grey display brightness of each frame relative to the mean of the frames
    just before and after it (NaN without a neighbour within JUMP_NEIGHBOUR_H).

    Earth's disc brightness changes smoothly over an hour; a jump that every
    colour shares and that the next frame undoes is the image processing."""
    g = df.loc[good].sort_values("t")
    grey = (g["D_red"] + g["D_green"] + g["D_blue"]) / 3.0
    dt_prev = g["t"].diff().dt.total_seconds() / 3600.0
    dt_next = -g["t"].diff(-1).dt.total_seconds() / 3600.0
    jump = grey / (0.5 * (grey.shift(1) + grey.shift(-1))) - 1.0
    jump[(dt_prev > JUMP_NEIGHBOUR_H) | (dt_next > JUMP_NEIGHBOUR_H)] = np.nan
    return jump.reindex(df.index)


def load(kind: str = "jpg") -> pd.DataFrame:
    """The light-curve table for one image size (``jpg`` or ``thumbs``)."""
    meta = epicdata.load_metadata()
    tab = pd.read_parquet(ph.table_paths(kind)[0])
    df = tab.merge(meta, on="image", how="inner").sort_values("t").reset_index(drop=True)
    size = 1080 if kind == "jpg" else 120
    df = add_geometry(df, size)
    df = add_fluxes(df)
    return quality(df)


def window(df: pd.DataFrame, start: str, end: str, *, ok_only: bool = True) -> pd.DataFrame:
    w = epicdata.between(df, start, end)
    return w[w["ok"]] if ok_only else w


def days_since(t: pd.Series, t0=None) -> np.ndarray:
    t0 = t.iloc[0] if t0 is None else pd.Timestamp(t0)
    return (t - t0).dt.total_seconds().to_numpy() / 86400.0
