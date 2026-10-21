"""The EPIC image catalogue: one row per natural-colour image.

Built from the per-day API files in ``data/raw/api/natural/`` and cached as
``data/cache/metadata.parquet``. Positions are geocentric J2000 vectors in km,
as the API gives them.

Two times describe each image. ``t_start`` (the API's ``date``) is when the
ten-band sequence began; the image name carries the time of EPIC's reference
band, almost always 288-289 s later. The natural-colour picture is built from
the 443, 551 and 680 nm bands taken inside that window, so we use the
reference time: the difference is far below anything the analysis resolves.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from config import (
    API_NATURAL_DIR, CACHE_DIR, EPIC_ARCHIVE, IMAGE_KINDS, NATURAL_DIR,
)

METADATA_PATH = CACHE_DIR / "metadata.parquet"


def _rows(path):
    for im in json.loads(path.read_text()):
        c = im["coords"]
        d, m, s = (c["dscovr_j2000_position"], c["lunar_j2000_position"],
                   c["sun_j2000_position"])
        q = c["attitude_quaternions"]
        yield {
            "image": im["image"], "version": im.get("version", ""),
            "t_start": im["date"],
            "lat": c["centroid_coordinates"]["lat"], "lon": c["centroid_coordinates"]["lon"],
            "dsc_x": d["x"], "dsc_y": d["y"], "dsc_z": d["z"],
            "moon_x": m["x"], "moon_y": m["y"], "moon_z": m["z"],
            "sun_x": s["x"], "sun_y": s["y"], "sun_z": s["z"],
            "q0": q["q0"], "q1": q["q1"], "q2": q["q2"], "q3": q["q3"],
        }


def build_metadata() -> pd.DataFrame:
    rows = [r for p in sorted(API_NATURAL_DIR.glob("*.json")) for r in _rows(p)]
    df = pd.DataFrame(rows)
    df["t_start"] = pd.to_datetime(df["t_start"])
    df["t"] = pd.to_datetime(df["image"].str[8:22], format="%Y%m%d%H%M%S")
    # A handful of images are listed twice (identical name): keep one.
    df = df.drop_duplicates("image").sort_values("t").reset_index(drop=True)
    df.to_parquet(METADATA_PATH)
    return df


def load_metadata(*, rebuild: bool = False) -> pd.DataFrame:
    if rebuild or not METADATA_PATH.exists():
        return build_metadata()
    newest = max(p.stat().st_mtime for p in API_NATURAL_DIR.glob("*.json"))
    if newest > METADATA_PATH.stat().st_mtime:
        return build_metadata()
    return pd.read_parquet(METADATA_PATH)


def between(df: pd.DataFrame, start: str, end: str) -> pd.DataFrame:
    """Images from the start of ``start`` to the end of ``end`` (dates, inclusive)."""
    t0 = pd.Timestamp(start)
    t1 = pd.Timestamp(end) + pd.Timedelta(days=1)
    return df[(df["t"] >= t0) & (df["t"] < t1)]


def image_url(kind: str, image: str) -> str:
    y, m, d = image[8:12], image[12:14], image[14:16]
    return f"{EPIC_ARCHIVE}/natural/{y}/{m}/{d}/{kind}/{image}.{IMAGE_KINDS[kind]}"


def image_path(kind: str, image: str):
    return NATURAL_DIR / kind / image[8:12] / image[12:14] / f"{image}.{IMAGE_KINDS[kind]}"


def unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def vectors(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Geocentric J2000 positions (km) of DSCOVR, the Moon and the Sun."""
    dsc = df[["dsc_x", "dsc_y", "dsc_z"]].to_numpy(float)
    moon = df[["moon_x", "moon_y", "moon_z"]].to_numpy(float)
    sun = df[["sun_x", "sun_y", "sun_z"]].to_numpy(float)
    return dsc, moon, sun
