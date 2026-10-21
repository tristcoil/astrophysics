"""Turning each picture of Earth into one point of light.

A natural-colour EPIC image is 8-bit RGB with space masked to exactly 0. For
every image we measure the disc (pixel count, centre, radius) and keep, per
colour channel, a 256-bin histogram of the disc's pixel values: once over the
whole sunlit disc and once over the inner disc (radius < ``INNER_FRACTION`` of
the disc radius, which drops the limb pixels where the Sun is within 12 deg of
the horizon, as Yang et al. 2018 did).

JPEG pixels are integers, so the histograms hold everything any per-pixel
transform needs: the disc mean of f(V) is ``hist @ f(v) / n``. That lets us
choose the display-to-linear transform afterwards (see ``transform.py``)
without re-reading ten thousand images.

The channels are EPIC's 680 nm (red), 551 nm (green) and 443 nm (blue) bands.
"""

from __future__ import annotations

import multiprocessing as mp
from pathlib import Path

import numpy as np
from PIL import Image

from config import CACHE_DIR
from transform import GAMMA_GRID

#: Inner-disc cut: cos(78 deg) = 0.208, and at the small phase angles EPIC
#: sees, the solar zenith angle is close to the viewing zenith angle, which
#: reaches 78 deg at sin(78 deg) = 0.978 of the disc radius.
INNER_FRACTION = float(np.sin(np.radians(78.0)))
#: A pixel belongs to the disc if its three channels sum to more than this
#: (out of 765). Space is 0; JPEG ringing beside the limb reaches a few units;
#: the sunlit limb itself is above 60 in every channel.
DISC_THRESHOLD = 24
CHANNELS = ("red", "green", "blue")
CHANNEL_NM = {"red": 680, "green": 551, "blue": 443}


def load_rgb(path: Path) -> np.ndarray:
    with Image.open(path) as im:
        return np.asarray(im.convert("RGB"))


def disc_mask(rgb: np.ndarray) -> np.ndarray:
    return rgb.astype(np.int16).sum(axis=-1) > DISC_THRESHOLD


def measure(rgb: np.ndarray) -> dict:
    """Disc geometry and per-channel histograms of one image."""
    ny, nx, _ = rgb.shape
    mask = disc_mask(rgb)
    n = int(mask.sum())
    yy, xx = np.nonzero(mask)
    cy, cx = float(yy.mean()), float(xx.mean())
    # Radius from the extent along both axes (robust to a thin unlit crescent
    # on one side, which shifts the centroid but not the opposite limb much).
    r_area = float(np.sqrt(n / np.pi))
    r_extent = 0.25 * float((xx.max() - xx.min() + 1) + (yy.max() - yy.min() + 1))
    yg, xg = np.ogrid[:ny, :nx]
    rr = np.hypot(xg - 0.5 * (xx.max() + xx.min()), yg - 0.5 * (yy.max() + yy.min()))
    inner = mask & (rr < INNER_FRACTION * r_extent)
    hist = np.zeros((2, 3, 256), dtype=np.uint32)
    for c in range(3):
        ch = rgb[..., c]
        hist[0, c] = np.bincount(ch[mask], minlength=256)
        hist[1, c] = np.bincount(ch[inner], minlength=256)
    out = {
        "nx": nx, "n_disc": n, "n_inner": int(inner.sum()),
        "cx": cx, "cy": cy, "r_area": r_area, "r_extent": r_extent,
        "x_min": int(xx.min()), "x_max": int(xx.max()), "y_min": int(yy.min()), "y_max": int(yy.max()),
        "n_sat": [int((rgb[..., c][mask] >= 254).sum()) for c in range(3)],
        "hist": hist,
    }
    # Disc means of (V/255)^gamma for each stored exponent, both regions.
    v = np.arange(256) / 255.0
    for g in GAMMA_GRID:
        f = v ** g
        for region, h in (("disc", hist[0]), ("inner", hist[1])):
            means = (h @ f) / np.maximum(h.sum(axis=-1), 1)
            for name, m in zip(CHANNELS, means):
                out[moment_column(region, name, g)] = float(m)
    return out


def moment_column(region: str, channel: str, gamma: float) -> str:
    """Column name of a stored disc mean, e.g. 'disc_red_g2.2'."""
    return f"{region}_{channel}_g{gamma:.1f}"


def _measure_path(args):
    image, path, keep_hist = args
    try:
        m = measure(load_rgb(path))
    except Exception as err:  # unreadable file: report, don't crash the pool
        return image, None, repr(err)
    if not keep_hist:
        m["hist"] = None
    return image, m, None


def _init_worker():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def measure_many(items: list[tuple[str, Path]], *, workers: int = 4, keep_hist: bool = True,
                 quiet: bool = False):
    """Measure many images in parallel. Returns (table rows, histograms, errors);
    histograms are (n, 2, 3, 256), or empty when ``keep_hist`` is False."""
    rows, hists, errors = [], [], {}
    ctx = mp.get_context("fork")  # Python 3.14 defaults to forkserver
    with ctx.Pool(workers, initializer=_init_worker) as pool:
        args = [(image, path, keep_hist) for image, path in items]
        for i, (image, m, err) in enumerate(pool.imap(_measure_path, args, chunksize=16)):
            if err is not None:
                errors[image] = err
                continue
            h = m.pop("hist")
            if keep_hist:
                hists.append(h)
            sat = m.pop("n_sat")
            m.update({f"n_sat_{c}": s for c, s in zip(CHANNELS, sat)})
            m["image"] = image
            rows.append(m)
            if not quiet and (i + 1) % 2000 == 0:
                print(f"  measured {i + 1}/{len(items)}", flush=True)
    empty = np.zeros((0, 2, 3, 256), np.uint32)
    return rows, (np.stack(hists) if hists else empty), errors


# ---------------------------------------------------------------------------
# Cached photometry tables
# ---------------------------------------------------------------------------

def table_paths(kind: str) -> tuple[Path, Path]:
    return CACHE_DIR / f"photometry_{kind}.parquet", CACHE_DIR / f"histograms_{kind}.npz"


def photometry_table(kind: str, images: list[str], *, workers: int = 4, quiet: bool = False):
    """Measure every listed image of one size that is on disk, reusing the cache.

    Returns (table, histograms or None, image order of the histograms)."""
    import pandas as pd

    import epicdata

    tab_path, hist_path = table_paths(kind)
    keep_hist = kind != "thumbs"  # 55,000 thumbnails: moments only
    old = pd.read_parquet(tab_path) if tab_path.exists() else None
    old_h, old_names = None, []
    if keep_hist and hist_path.exists():
        z = np.load(hist_path, allow_pickle=False)
        old_h, old_names = z["hist"], list(z["image"])
    done = set(old["image"]) if old is not None else set()
    todo = [(im, epicdata.image_path(kind, im)) for im in images
            if im not in done and epicdata.image_path(kind, im).exists()]
    if todo:
        if not quiet:
            print(f"photometry ({kind}): measuring {len(todo)} images", flush=True)
        rows, hists, errors = measure_many(todo, workers=workers, keep_hist=keep_hist, quiet=quiet)
        if errors and not quiet:
            print(f"  {len(errors)} unreadable: {list(errors)[:3]}")
        new = pd.DataFrame(rows)
        tab = new if old is None else pd.concat([old, new], ignore_index=True)
        tab = tab.drop_duplicates("image").reset_index(drop=True)
        tab.to_parquet(tab_path)
        if keep_hist:
            names = old_names + [r["image"] for r in rows]
            allh = hists if old_h is None else np.concatenate([old_h, hists])
            np.savez_compressed(hist_path, hist=allh, image=np.array(names))
            old_h, old_names = allh, names
    else:
        tab = old
    wanted = set(images)
    tab = tab[tab["image"].isin(wanted)].reset_index(drop=True) if tab is not None else None
    return tab, old_h, old_names


def disc_means(hist: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Disc mean of f(V) for each image and channel.

    hist: (..., 256) counts; f: (256,) or (3, 256) values of the transform per
    pixel value (per channel if 2-D). Returns hist's shape without the last axis.
    """
    n = hist.sum(axis=-1)
    if f.ndim == 1:
        s = hist @ f
    else:
        s = np.einsum("...cv,cv->...c", hist, f)
    return s / np.maximum(n, 1)
