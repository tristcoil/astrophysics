"""The Sun as a star: calibrating the scaling relations with our own code.

The scaling relations are anchored on the Sun's numax and Dnu, and every
pipeline measures the Sun's numax a little differently (Pinsonneault et al.
2025 list effective solar values of 3057-3090 uHz). A mass goes as numax^3,
so a 1% mismatch between how the stars and the Sun were measured is a 3%
error in mass. We therefore run the *same* code on 25 years of SOHO/VIRGO
photometry of the whole Sun, one calendar year at a time, and use the result
as our solar reference.

VIRGO's green sun photometer (500 nm) samples the Sun's brightness every
minute; Kepler's band is broader (~420-900 nm). Mode amplitudes depend on
wavelength, but numax and Dnu do not.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import lightcurves as L
from fetch_data import VIRGO_GREEN_PATH
from stars import global_parameters

SUN_RESULTS = C.CACHE_DIR / "sun.json"
SUN_ARRAYS = C.CACHE_DIR / "sun.npz"
#: Calendar years analysed. 1998 lost half its data (SOHO was out of contact
#: from June to October), and 2023 ends in April; both are kept, flagged by
#: their duty cycle.
YEARS = tuple(range(1996, 2023))
#: Settings for a Sun-like star: the fit starts at 30 uHz (well above the
#: 1-day detrending filter) and numax is searched from 1000 to 5000 uHz. A
#: fully blind search (300-8000 uHz) locked onto an instrumental hump near
#: 5.5 mHz in 2016 (Huber et al. 2009 note a ~6 mHz artefact in VIRGO
#: photometry), so for the Sun the search is not blind; it does not need to be.
SUN_SEARCH = (1000.0, 5000.0, 24)
SUN_FIT_MIN = 30.0


def virgo_series() -> tuple[np.ndarray, np.ndarray, float]:
    """VIRGO green channel: time (days since 1996-01-23 TAI) and flux (ppm)."""
    from astropy.io import fits

    d = fits.getdata(VIRGO_GREEN_PATH)
    t = d[:, 0] / 86400.0
    f = d[:, 1]
    return t, f, 1996.0 + 22.0 / 365.25  # 1996-01-23 as a decimal year


def year_segment(t, f, year: int) -> L.TimeSeries:
    """One calendar year of VIRGO, prepared like the Kepler data.

    VIRGO's level-2 series is already in ppm around zero with gaps as NaN;
    we add 1 to make it a relative flux, cut the year into day-long-gap
    separated pieces and run the same detrend-and-clip recipe (1-day trend).
    """
    t0 = (np.datetime64(f"{year}-01-01") - np.datetime64("1996-01-23")).astype(float)
    t1 = (np.datetime64(f"{year + 1}-01-01") - np.datetime64("1996-01-23")).astype(float)
    m = (t >= t0) & (t < t1) & np.isfinite(f)
    tt, ff = t[m], 1.0 + f[m] * 1e-6
    # cadence index on the regular 60-s grid (the one 90-s step in March 1996
    # shifts later samples by half a cadence; floor(x + 0.25) absorbs it)
    cad = np.floor((tt - t0) * 86400.0 / C.VIRGO_CADENCE_S + 0.25).astype(np.int64)
    return L.prepare_segments([(tt, ff, cad)], cadence_s=C.VIRGO_CADENCE_S,
                              trend_sigma_days=C.TREND_SIGMA_DAYS_SC, label=f"Sun {year}")


_VIRGO = {}


def _year_worker(y):
    from threadpoolctl import threadpool_limits

    t, f = _VIRGO["t"], _VIRGO["f"]
    with threadpool_limits(1):
        ts = year_segment(t, f, y)
        out, arr = global_parameters(ts, seed=y, search=SUN_SEARCH, fit_nu_min=SUN_FIT_MIN,
                                     period_spacing=False, eps_width=0.5, n_perturb=100,
                                     dnu_method="acf")
    return y, out, arr


def analyse_sun(*, force: bool = False, quiet: bool = False, workers: int = 6) -> dict:
    """Every calendar year of VIRGO through the giants' pipeline (cached)."""
    import multiprocessing as mp

    if SUN_RESULTS.exists() and not force:
        return json.loads(SUN_RESULTS.read_text())
    t, f, _ = virgo_series()
    _VIRGO["t"], _VIRGO["f"] = t, f  # inherited by the forked workers
    years, arrays = {}, {}
    with mp.get_context("fork").Pool(workers) as pool:
        for y, out, arr in pool.imap_unordered(_year_worker, YEARS):
            years[y] = out
            if not quiet:
                print(f"    Sun {y}: duty {out['duty_cycle']:.2f}  numax {out['numax']:.0f} "
                      f"+- {out['numax_err']:.0f}  Dnu {out['dnu']:.3f} +- {out['dnu_err']:.3f} uHz",
                      flush=True)
            if y in (2008, 2014):  # solar minimum and maximum, for figures
                for k, v in arr.items():
                    arrays[f"y{y}_{k}"] = v
    years = {y: years[y] for y in sorted(years)}
    good = [y for y in YEARS if years[y]["duty_cycle"] > 0.8]
    nm = np.array([years[y]["numax"] for y in good])
    dn = np.array([years[y]["dnu"] for y in good])
    res = {
        "years": years,
        "years_used": good,
        "numax_sun": float(np.mean(nm)),
        "numax_sun_scatter": float(np.std(nm, ddof=1)),
        "numax_sun_err": float(np.std(nm, ddof=1) / np.sqrt(len(nm))),
        "dnu_sun": float(np.mean(dn)),
        "dnu_sun_scatter": float(np.std(dn, ddof=1)),
        "dnu_sun_err": float(np.std(dn, ddof=1) / np.sqrt(len(dn))),
    }
    np.savez_compressed(SUN_ARRAYS, **arrays)
    SUN_RESULTS.write_text(json.dumps(res, indent=1, default=float))
    return res
