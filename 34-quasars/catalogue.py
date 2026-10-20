"""The DR16Q quasar catalogue joined to Wu & Shen's measurements, with derived columns.

Adds, for each of the 750,414 quasars:
- ``programme``: which SDSS survey programme put a fibre on it (from DR16Q's
  target bits), the root of the catalogue's uneven selection;
- dereddened ugriz PSF magnitudes and colours;
- M_i(z=2), the luminosity measure of Richards et al. (2006), in D3's cosmology;
- black hole masses, bolometric luminosities and Eddington ratios (masses.py).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import masses
import physics
from config import DR16Q_TABLE, WUSHEN_TABLE

#: Programme labels, in the priority used when a quasar carries several flags.
PROGRAMMES = ("eBOSS CORE", "BOSS CORE", "SDSS-I/II", "Other BOSS/eBOSS")
PROGRAMME_NOTES = {
    "eBOSS CORE": "uniform XDQSOz + WISE selection to g or r < 22, aimed at 0.9 < z < 2.2 (2014-2019)",
    "BOSS CORE": "uniform XDQSO selection to g < 22.0 or r < 21.85, aimed at z > 2.15 (2009-2014)",
    "SDSS-I/II": "ugriz colour + FIRST radio selection to i < 19.1 (z < 3) or 20.2 (2000-2008)",
    "Other BOSS/eBOSS": "BONUS, variability, known-quasar re-observations and ancillary programmes",
}

#: eBOSS_TARGET1 bit 10 = QSO1_EBOSS_CORE; BOSS_TARGET1 bit 40 = QSO_CORE_MAIN,
#: bit 10 = QSO_CORE (the DR9-era name of the same uniform sample).
_EBOSS_CORE = np.uint64(1 << 10)
_BOSS_CORE = np.uint64((1 << 40) | (1 << 10))
#: The last SDSS-I/II spectra predate MJD 55000 (June 2009); BOSS began that December.
SDSS_LAST_MJD = 55000


#: the Wu & Shen columns the analysis uses (the cached table keeps more)
USED_COLUMNS = (["SDSS_NAME", "PLATE", "MJD", "FIBERID", "RA", "DEC", "IF_BOSS_SDSS", "Z_DR16Q",
                 "Z_FIT", "Z_SYS", "Z_SYS_ERR", "EBV", "SN_MEDIAN_ALL", "LOGL1350", "LOGL3000",
                 "LOGL5100", "HBETA_BR_2", "HBETA_BR_4", "MGII_BR_2", "MGII_BR_4", "CIV_1", "CIV_2",
                 "CIV_4", "LOGMBH_HB_ERR", "LOGMBH_MGII_ERR", "LOGMBH_CIV_ERR", "LOGMBH", "LOGLBOL",
                 "LOGLEDD_RATIO", "ZSYS_LINES_4", "ZSYS_LINES_ERR_4", "FHOST_5100", "BAL_PROB", "M_I",
                 "FIRST_MATCHED", "HBETA_BR_ERR_4", "MGII_BR_ERR_4", "CIV_ERR_4", "LOGL5100_ERR",
                 "LOGL3000_ERR", "LOGL1350_ERR"]
                + [f"PSFMAG_{j}" for j in range(5)] + [f"EXTINCTION_{j}" for j in range(5)])


def programme_labels(q: pd.DataFrame) -> np.ndarray:
    e1 = q["eBOSST1"].to_numpy(np.uint64)
    b1 = q["Btarg1"].to_numpy(np.uint64)
    lab = np.full(len(q), PROGRAMMES[3], dtype=object)
    sdss = q["MJD"].to_numpy() < SDSS_LAST_MJD
    lab[sdss] = PROGRAMMES[2]
    lab[(b1 & _BOSS_CORE) > 0] = PROGRAMMES[1]
    lab[(e1 & _EBOSS_CORE) > 0] = PROGRAMMES[0]
    return lab


def load(kcorr=None) -> pd.DataFrame:
    """The joined catalogue. ``kcorr`` maps z -> K(z) for M_i(z=2) (Richards et al. by default)."""
    ws = pd.read_parquet(WUSHEN_TABLE, columns=USED_COLUMNS)
    q = pd.read_parquet(DR16Q_TABLE)
    q = q.rename(columns={"Plate": "PLATE", "MJD": "MJD", "Fiber": "FIBERID"})
    keep = ["PLATE", "MJD", "FIBERID", "z", "r_z", "Btarg1", "eBOSST0", "eBOSST1", "eBOSST2",
            "Atarg1", "Atarg2"]
    cat = ws.merge(q[keep], on=["PLATE", "MJD", "FIBERID"], how="left", validate="one_to_one",
                   indicator=True)
    if (cat["_merge"] != "both").any():
        raise RuntimeError(f"{int((cat['_merge'] != 'both').sum())} Wu & Shen rows not in DR16Q")
    cat = cat.drop(columns="_merge")
    cat["programme"] = programme_labels(cat.rename(columns={}))

    # Redshift: Wu & Shen's systemic redshift where valid, else DR16Q's.
    zsys = cat["Z_SYS"].to_numpy()
    cat["zbest"] = np.where(zsys > 0, zsys, cat["Z_DR16Q"].to_numpy())

    # Dereddened PSF magnitudes (asinh mags; fine for the bright objects we use them for).
    for j, b in enumerate("ugriz"):
        cat[f"{b}0"] = cat[f"PSFMAG_{j}"] - cat[f"EXTINCTION_{j}"]
    for a, b in zip("ugri", "griz"):
        cat[f"{a}-{b}"] = cat[f"{a}0"] - cat[f"{b}0"]

    if kcorr is not None:
        z = cat["zbest"].to_numpy()
        ok = z > 0
        mi = np.full(len(cat), np.nan)
        mi[ok] = physics.absolute_mag(cat.loc[ok, "i0"].to_numpy(), z[ok], kcorr, "d3")
        cat["Mi_z2"] = mi
        cat["logL2500"] = physics.luminosity_2500(mi)

    m = masses.recompute(cat)
    for c in m.columns:
        cat[c] = m[c]
    return cat
