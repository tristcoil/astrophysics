"""GALAH DR4: load the catalogue once, keep what we use, define samples.

The 723 MB catalogue is read column by column (memory-mapped) and the stars
passing the survey's own recommendations are cached as a compact .npz:

* field programmes only (config.FIELD_SURVEYS),
* ``flag_sp == 0`` (no warning on the stellar parameters),
* ``flag_fe_h == 0`` and S/N per pixel in the red camera > 30.

An abundance [X/Fe] is used only where ``flag_x_fe == 0`` (a detection, not
an upper limit, with no fitting warning); otherwise it is NaN in the cache.
"""

from __future__ import annotations

import numpy as np

from config import (CACHE_DIR, DWARFS, ELEMENTS, EU_DIR, FIELD_SURVEYS, GALAH_DIR,
                    GALAH_FILES, GIANTS, GSE_LZ, GSE_SQRT_JR, SNR_MIN)

CACHE = CACHE_DIR / "galah_dr4_clean.npz"

_BASIC = ["sobject_id", "gaiadr3_source_id", "ra", "dec", "teff", "logg", "fe_h", "e_fe_h",
          "snr_px_ccd2", "snr_px_ccd3", "vsini", "parallax", "parallax_error", "ruwe",
          "rv_comp_1", "age", "mass", "phot_g_mean_mag", "bp_rp"]
_DYN = ["L_Z", "J_R", "J_Z", "Energy", "ecc", "zmax", "R_Rzphi", "z_Rzphi", "vR_Rzphi",
        "vT_Rzphi", "vz_Rzphi", "R_peri", "R_ap"]


def _build_cache() -> None:
    from astropy.io import fits

    path = GALAH_DIR / GALAH_FILES["allstar"]
    with fits.open(path, memmap=True) as h:
        d = h[1].data
        n_all = len(d)
        survey = np.char.strip(np.asarray(d["survey_name"]).astype(str))
        keep = (np.isin(survey, FIELD_SURVEYS) & (np.asarray(d["flag_sp"]) == 0)
                & (np.asarray(d["flag_fe_h"]) == 0) & (np.asarray(d["snr_px_ccd3"]) > SNR_MIN))
        counts = dict(n_catalogue=int(n_all), n_field=int(np.isin(survey, FIELD_SURVEYS).sum()),
                      n_clean=int(keep.sum()))
        out = {c: np.asarray(d[c])[keep] for c in _BASIC}
        for el in ELEMENTS:
            if el == "Fe":
                continue
            e = el.lower()
            ok = np.asarray(d[f"flag_{e}_fe"])[keep] == 0
            out[f"{e}_fe"] = np.where(ok, np.asarray(d[f"{e}_fe"])[keep], np.nan).astype(np.float32)
            out[f"e_{e}_fe"] = np.asarray(d[f"e_{e}_fe"])[keep].astype(np.float32)
            # DR4's own (unreliable) europium is kept for the comparison figure.
        out["flag_eu_fe"] = np.asarray(d["flag_eu_fe"])[keep]
        out["eu_fe_dr4_raw"] = np.asarray(d["eu_fe"])[keep].astype(np.float32)

    # Orbits (galah_dr4_vac_dynamics, same row order by sobject_id).
    with fits.open(GALAH_DIR / GALAH_FILES["dynamics"], memmap=True) as h:
        dd = h[1].data
        sid = np.asarray(dd["sobject_id"])
        order = np.argsort(sid)
        pos = np.searchsorted(sid[order], out["sobject_id"])
        pos = np.clip(pos, 0, sid.size - 1)
        found = sid[order][pos] == out["sobject_id"]
        for c in _DYN:
            v = np.asarray(dd[c])[order][pos].astype(np.float64)
            out[c] = np.where(found, v, np.nan)
    counts["n_with_orbits"] = int(np.isfinite(out["L_Z"]).sum())
    for k, v in list(out.items()):
        if v.dtype.kind == "f" and k not in ("ra", "dec") and k not in _DYN:
            out[k] = v.astype(np.float32)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, **out, **{f"_count_{k}": v for k, v in counts.items()})


def load() -> dict:
    """The clean GALAH sample as a dict of arrays (+ '_counts')."""
    if not CACHE.exists():
        _build_cache()
    with np.load(CACHE) as z:
        d = {k: z[k] for k in z.files if not k.startswith("_count_")}
        d["_counts"] = {k[7:]: int(z[k]) for k in z.files if k.startswith("_count_")}
    d["mg_h"] = d["mg_fe"] + d["fe_h"]
    d["sqrt_jr"] = np.sqrt(np.clip(d["J_R"], 0, None))
    return d


def in_range(x, lo_hi):
    return (x > lo_hi[0]) & (x < lo_hi[1])


def dwarfs(d):
    return in_range(d["logg"], DWARFS["logg"]) & in_range(d["teff"], DWARFS["teff"])


def giants(d):
    return in_range(d["logg"], GIANTS["logg"]) & in_range(d["teff"], GIANTS["teff"])


def gse_dynamical(d):
    """Feuillet et al. (2021) / Buder et al. (2022): -500 < L_Z < 500 kpc km/s,
    30 < sqrt(J_R) < 55 sqrt(kpc km/s). Applied to GALAH's own orbits."""
    return in_range(d["L_Z"], GSE_LZ) & in_range(d["sqrt_jr"], GSE_SQRT_JR)


def disc_kinematics(d):
    """Stars on disc-like orbits: prograde, not eccentric and not far from the plane."""
    return (d["L_Z"] > 0) & (d["ecc"] < 0.5) & (d["zmax"] < 3.0)


# ---------------------------------------------------------------------------
# Europium (Kane et al. 2026 value-added catalogue)
# ---------------------------------------------------------------------------


def load_europium(d) -> dict:
    """Attach [Eu/H] from the europium VAC to the GALAH sample, by sobject_id.

    * ``eu_cnn``: CNN predictions with eu_flag == 0 (all Kane et al.'s
      recommended cuts: 4000 < Teff < 6600 K, [Fe/H] > -1.7, clean continuum,
      a visible line);
    * ``eu_golden``: their "golden sample" (giants with eu_flag == 0, plus
      Korg syntheses for -1.7 < [Fe/H] < -1);
    * ``eu_dr3``: GALAH DR3's own [Eu/H] for the 52,147 stars used to train
      the network (the labels the CNN learnt to reproduce).
    All are returned as [Eu/Fe] = [Eu/H] - [Fe/H]_DR4.
    """
    from astropy.io import fits

    def attach(fname, col, flag=None, extra=None):
        with fits.open(EU_DIR / fname, memmap=True) as h:
            t = h[1].data
            sid = np.asarray(t["sobject_id"])
            val = np.asarray(t[col], float)
            ok = np.isfinite(val)
            if flag is not None:
                ok &= np.asarray(t[flag]) == 0
            if extra is not None:
                ok &= extra(t)
            sid, val = sid[ok], val[ok]
        order = np.argsort(sid)
        sid, val = sid[order], val[order]
        pos = np.clip(np.searchsorted(sid, d["sobject_id"]), 0, max(sid.size - 1, 0))
        hit = sid[pos] == d["sobject_id"]
        return np.where(hit, val[pos] - d["fe_h"], np.nan).astype(np.float32)

    eu = {
        "eu_cnn": attach("cnn_predictions.fits", "eu_h_pred", flag="eu_flag"),
        "eu_golden": attach("golden_sample.fits", "eu_h"),
        "eu_korg": attach("korg_secondary_catalogue.fits", "eu_h_korg", flag="eu_flag",
                          extra=lambda t: np.asarray(t["lower_limit_flag"]) == 0),
        "eu_dr3": attach("training_and_validation_dataset.fits", "eu_h_DR3"),
    }
    # The golden sample, as Kane et al. define it: CNN giants at [Fe/H] >= -1,
    # Korg below. Our "best" europium uses it where available, then DR3's own.
    best = np.where(np.isfinite(eu["eu_golden"]), eu["eu_golden"], eu["eu_dr3"])
    eu["eu_best"] = best.astype(np.float32)
    # Coverage: stars the CNN catalogue looked at and whose spectra pass its
    # non-abundance checks (stellar parameters, continuum). Among these, the
    # line-depth check (the only abundance-dependent one) passes for >88% of
    # giants above [Fe/H] = -1, so medians are not biased by non-detections.
    eu["cov_cnn"] = np.isfinite(attach("cnn_predictions.fits", "eu_h_pred",
                                       extra=lambda t: (np.asarray(t["stellar_params_flag"]) == 0)
                                       & (np.asarray(t["continuum_flux_flag"]) == 0)))
    return eu
