"""El-Badry, Rix & Heintz (2021) catalogue: slim copies and published samples.

The Zenodo file is a 1.42 GB gzipped FITS table (1,817,594 pairs x 217
columns, ~2.5 GB uncompressed). Reading it straight from the .gz would hold
all of it in memory, so we decompress once to a scratch file under
data/derived, memory-map it, keep two slim copies and delete the scratch file:

* ``elbadry_all``: every pair, a few columns (sky, distance, separation,
  chance-alignment ratio R, binary type), for the overview figures;
* ``elbadry_near``: pairs with parallax1 > 3 mas (closer than ~330 pc) and
  every column the kinematic samples need, joined row by row to Hwang, Ting
  & Zakamska's (2022) eccentricity estimates (same row order, checked by
  source_id).

The shifted chance-alignment catalogue (517,993 pairs) is slimmed the same
way.
"""

from __future__ import annotations

import gzip
import shutil

import numpy as np
from astropy.io import fits
from astropy.table import Table

from config import DERIVED_DIR, RAW_DIR

STAR_COLS = ["source_id", "ra", "dec", "l", "b", "parallax", "parallax_error", "pmra", "pmra_error",
             "pmdec", "pmdec_error", "ra_dec_corr", "ra_parallax_corr", "ra_pmra_corr", "ra_pmdec_corr",
             "dec_parallax_corr", "dec_pmra_corr", "dec_pmdec_corr", "parallax_pmra_corr",
             "parallax_pmdec_corr", "pmra_pmdec_corr", "phot_g_mean_mag", "phot_g_mean_flux_over_error",
             "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp", "phot_bp_rp_excess_factor", "ruwe",
             "ipd_frac_multi_peak", "ipd_gof_harmonic_amplitude", "astrometric_chi2_al",
             "astrometric_n_good_obs_al", "astrometric_excess_noise", "astrometric_params_solved",
             "visibility_periods_used", "dr2_radial_velocity", "dr2_radial_velocity_error"]
PAIR_COLS = ["pairdistance", "sep_AU", "binary_type", "Sigma18", "R_chance_align"]
ALL_COLS = ["source_id1", "source_id2", "ra1", "dec1", "l1", "b1", "parallax1", "parallax2",
            "parallax_error1", "parallax_error2", "phot_g_mean_mag1", "phot_g_mean_mag2", "bp_rp1",
            "bp_rp2", "pmra1", "pmdec1", "sep_AU", "pairdistance", "R_chance_align", "binary_type"]


def _decompressed(name: str):
    src = RAW_DIR / name
    tmp = DERIVED_DIR / (name.replace(".gz", "") + ".tmp")
    if not tmp.exists():
        print(f"  decompressing {name} to a scratch file ...", flush=True)
        with gzip.open(src, "rb") as fin, open(tmp.with_suffix(".part"), "wb") as fout:
            shutil.copyfileobj(fin, fout, 2**24)
        tmp.with_suffix(".part").rename(tmp)
    return tmp


def _slim(name: str, out_all: str, out_near: str, near_plx: float = 3.0, with_hwang: bool = False):
    tmp = _decompressed(name)
    with fits.open(tmp, memmap=True) as h:
        d = h[1].data
        n = len(d)
        cols = [c.name for c in h[1].columns]
        t_all = Table({c: np.asarray(d[c]) for c in ALL_COLS if c in cols})
        t_all.write(DERIVED_DIR / f"{out_all}.fits", overwrite=True)
        near = np.asarray(d["parallax1"]) > near_plx
        keep = [f"{c}{w}" for c in STAR_COLS for w in (1, 2)] + PAIR_COLS
        t_near = Table({c: np.asarray(d[c][near]) for c in keep if c in cols})
        t_near["row"] = np.flatnonzero(near)
        if with_hwang:
            # read with astropy.io.fits: the file's header carries astropy-table
            # metadata naming columns it does not have, which Table.read rejects
            with fits.open(RAW_DIR / "hwang2022" / "hwang2022_eccentricities.fits", memmap=True) as hh:
                hw = hh[1].data
                assert len(hw) == n
                sel = {c: np.asarray(hw[c][near]) for c in
                       ("source_id1", "source_id2", "vr_angle", "vr_angle_error", "dpm_sig", "alpha",
                        "e", "e0", "e1")}
            if not (np.array_equal(np.asarray(sel["source_id1"]), np.asarray(t_near["source_id1"]))
                    and np.array_equal(np.asarray(sel["source_id2"]), np.asarray(t_near["source_id2"]))):
                raise RuntimeError("Hwang et al. rows do not line up with El-Badry et al.'s")
            for c in ("vr_angle", "vr_angle_error", "dpm_sig", "alpha", "e", "e0", "e1"):
                t_near[f"hw_{c}"] = np.asarray(sel[c])
        t_near.write(DERIVED_DIR / f"{out_near}.fits", overwrite=True)
    tmp.unlink()
    return n


def ensure_slim(quiet: bool = False):
    if not (DERIVED_DIR / "elbadry_near.fits").exists():
        _slim("all_columns_catalog.fits.gz", "elbadry_all", "elbadry_near", with_hwang=True)
    if not (DERIVED_DIR / "elbadry_shift_near.fits").exists():
        _slim("all_columns_catalog_shift.fits.gz", "elbadry_shift_all", "elbadry_shift_near")


def load(which: str = "near") -> Table:
    ensure_slim()
    return Table.read(DERIVED_DIR / f"elbadry_{which}.fits")


if __name__ == "__main__":
    import time
    t0 = time.time()
    ensure_slim()
    for w in ("all", "near", "shift_all", "shift_near"):
        t = Table.read(DERIVED_DIR / f"elbadry_{w}.fits")
        print(w, len(t), len(t.colnames), f"{time.time() - t0:.0f} s")
