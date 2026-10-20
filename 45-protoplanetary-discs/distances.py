"""Distances to the host stars from Gaia DR3 parallaxes.

Converting an angle on the sky to a size in au needs the distance: 1 arcsec at 100 pc is
100 au. DSHARP used Gaia DR2 (2018); this project uses DR3 (2022), corrected for the
parallax zero point of Lindegren et al. (2021, the ``gaiadr3_zeropoint`` package), with
the formal errors inflated (factor 1.3, times the RUWE when above 1, plus 10.3 micro-arcsec
in quadrature, after Maiz Apellaniz et al. 2021) because several hosts are binaries or
have discs that disturb the photocentre. A parallax worse than 5% is replaced by the mean
of the same cloud's DSHARP members (this applies to Elias 27 only).
"""

from __future__ import annotations

import csv
import warnings
from functools import lru_cache

import numpy as np

from config import (DISTANCE_MAX_FRAC_ERROR, GAIA_PATH, HLTAU, HLTAU_DISTANCE_PC, PARALLAX_ERROR_FLOOR_FACTOR,
                    PARALLAX_SYSTEMATIC_MAS, TARGETS)


# Median parallax offset of quasars in Gaia (E)DR3 (Lindegren et al. 2021), used where the
# zero-point calibration does not apply, with this much extra uncertainty.
GLOBAL_ZERO_POINT_MAS = -0.017
GLOBAL_ZERO_POINT_UNC_MAS = 0.015


def _f(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


@lru_cache(maxsize=None)
def gaia_distances() -> dict[str, dict]:
    from zero_point import zpt  # the gaiadr3_zeropoint package installs as "zero_point"

    zpt.load_tables()
    rows = {r["key"]: r for r in csv.DictReader(open(GAIA_PATH))}
    out = {}
    for t in TARGETS:
        r = rows[t.key]
        plx, err, ruwe = _f(r["parallax"]), _f(r["parallax_error"]), _f(r["ruwe"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # the package wants numpy arrays (Python scalars trip numpy's NEP 50 casting rules)
            arr = lambda v: np.array([v])  # noqa: E731
            zp = float(zpt.get_zpt(arr(_f(r["phot_g_mean_mag"])), arr(_f(r["nu_eff_used_in_astrometry"])),
                                   arr(_f(r["pseudocolour"])), arr(_f(r["ecl_lat"])),
                                   arr(int(r["astrometric_params_solved"])), _warnings=False)[0])
        zp_note = "Lindegren et al. 2021"
        extra = 0.0
        if not np.isfinite(zp):  # outside the calibration's colour/magnitude range
            zp, extra, zp_note = GLOBAL_ZERO_POINT_MAS, GLOBAL_ZERO_POINT_UNC_MAS, "global quasar offset"
        plx_c = plx - zp
        err_ext = float(np.sqrt((PARALLAX_ERROR_FLOOR_FACTOR * max(1.0, ruwe) * err) ** 2
                                + PARALLAX_SYSTEMATIC_MAS**2 + extra**2))
        d = 1000.0 / plx_c
        out[t.key] = dict(source_id=int(r["source_id"]), parallax=plx, parallax_error=err, ruwe=ruwe,
                          gmag=_f(r["phot_g_mean_mag"]), zero_point=zp, zero_point_source=zp_note,
                          parallax_corr=plx_c,
                          parallax_error_ext=err_ext, d_parallax=d, d_parallax_err=d * err_ext / plx_c,
                          d_bailer_jones=_f(r["r_med_geo"]), region=t.region)
    # the adopted distance
    for key, s in out.items():
        if s["parallax_error_ext"] / s["parallax_corr"] <= DISTANCE_MAX_FRAC_ERROR:
            s.update(distance=s["d_parallax"], distance_err=s["d_parallax_err"], distance_source="Gaia DR3 parallax")
    for key, s in out.items():
        if "distance" in s:
            continue
        mates = [m for k, m in out.items() if k != key and m["region"] == s["region"] and "distance" in m
                 and m["distance_source"].startswith("Gaia DR3 parallax")]
        w = np.array([1.0 / m["distance_err"] ** 2 for m in mates])
        dm = float(np.sum(w * [m["distance"] for m in mates]) / np.sum(w))
        spread = float(np.std([m["distance"] for m in mates], ddof=1))
        s.update(distance=dm, distance_err=spread,
                 distance_source=f"mean of {len(mates)} {s['region']} members (own parallax "
                                 f"{s['parallax']:.2f} +/- {s['parallax_error']:.2f} mas is too uncertain)")
    out[HLTAU.key] = dict(distance=HLTAU_DISTANCE_PC[0], distance_err=HLTAU_DISTANCE_PC[1],
                          distance_source="VLBA parallaxes of L1551 (Galli et al. 2018); Gaia DR3 has no parallax",
                          region=HLTAU.region)
    return out


def distance(key: str) -> float:
    return gaia_distances()[key]["distance"]


if __name__ == "__main__":
    import literature as L

    dr2 = L.dsharp1_stars()
    for k, s in gaia_distances().items():
        old = dr2[k].distance.value if k in dr2 else np.nan
        print(f"{k:9s} d = {s['distance']:6.1f} +/- {s['distance_err']:4.1f} pc  (DSHARP/DR2 {old:5.0f}; "
              f"ZP {s.get('zero_point', np.nan) * 1e3:+.1f} uas, RUWE {s.get('ruwe', np.nan):.2f})  {s['distance_source']}")
