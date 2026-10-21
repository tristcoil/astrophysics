"""Load the cached archives into analysis-ready velocity sets.

Four spectrographs, three decades:

* **ELODIE** (1.93-m, Haute-Provence, 1994-2006): the discovery instrument.
  Velocities come from the archive's cross-correlation of each spectrum with a
  K0-star mask; the archive quotes no per-spectrum uncertainty, so each gets a
  nominal 1 m/s and the real scatter is fitted as jitter.
* **Lick/Hamilton** (3-m and 0.6-m CAT, 1987-2011): iodine-cell velocities.
  Each CCD dewar gets its own zero point; everything before 1995, when the
  precision was ~10 m/s and four dewars alternated, is one more.
* **Keck/HIRES** (10-m, 1996-2014): iodine-cell velocities. The August 2004
  CCD upgrade splits it into two instruments.
* **HARPS** (3.6-m, La Silla, 2003-): SERVAL velocities corrected for nightly
  zero points (RVBank). The May 2015 fibre exchange splits it in two.

Times are put on one scale, BJD(TDB). Lick's are Julian dates of the
observation, converted here with the star's position and the observatory's;
ELODIE's are barycentric JDs on the UTC scale. Every cut is logged, so the
README can say exactly how many points each step removed.
"""

from __future__ import annotations

import gzip
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from config import (
    ELODIE_IMATYP_PREFIX, ELODIE_MASK, FLARE_HALPHA_NSIGMA, HARPS_FIBRE_BJD,
    HIRES_UPGRADE_BJD, KECK, KECK_TALOR, LICK, LICK_EARLY_END_BJD, LINE_SHAPE_NSIGMA,
    MAX_ERROR_OVER_MEDIAN, RVBANK_GOOD_FLAG, STARS, Star,
)
from fetch_data import elodie_ccf_path, elodie_times_path, parse_elodie_ccf, rvbank_path
from rv_fit import RVData

#: Nominal per-spectrum error given to ELODIE velocities (see module docstring).
ELODIE_NOMINAL_ERR = 1.0

#: Activity indicators carried along from HARPS RVBank (column -> our name).
HARPS_ACTIVITY = {
    "Halpha": "Halpha", "e_Halpha": "Halpha_err",
    "CRX": "CRX", "e_CRX": "CRX_err",
    "dLW": "dLW", "e_dLW": "dLW_err",
    "NaD1": "NaD1", "NaD2": "NaD2",
    "RHKp": "RHK", "dRHKp": "RHK_err",
    "FWHM_DRS": "FWHM", "BIS": "BIS", "CONTRAST_DRS": "contrast",
}


@dataclass
class CutLog:
    """Point counts after each cut, per source, for the README."""

    steps: list[tuple[str, str, int]] = field(default_factory=list)

    def add(self, source: str, label: str, n: int) -> None:
        self.steps.append((source, label, int(n)))

    def as_dict(self) -> dict:
        out: dict[str, list] = {}
        for source, label, n in self.steps:
            out.setdefault(source, []).append({"step": label, "n": n})
        return out


# ---------------------------------------------------------------------------
# Time scales
# ---------------------------------------------------------------------------

def jd_utc_to_bjd_tdb(jd, star: Star, site: str):
    """Geocentric/topocentric JD(UTC) -> barycentric JD(TDB) for a star.

    Adds the light-travel time across the Earth's orbit (up to +-8.3 min)
    and the UTC -> TDB offset (~1 min).
    """
    from astropy import units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        loc = EarthLocation.of_site(site)
        t = Time(np.asarray(jd, float), format="jd", scale="utc", location=loc)
        coord = SkyCoord(star.ra_deg * u.deg, star.dec_deg * u.deg)
        return (t.tdb + t.light_travel_time(coord, "barycentric")).jd


def bjd_utc_to_tdb(bjd):
    """A barycentric JD recorded on the UTC scale -> TDB (adds TT-UTC)."""
    from astropy.time import Time

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Time(np.asarray(bjd, float), format="jd", scale="utc").tdb.jd


# ---------------------------------------------------------------------------
# Per-instrument loaders -> tidy DataFrames with columns t, rv, err, inst
# ---------------------------------------------------------------------------

def _read_keck_table(version: str = "butler") -> pd.DataFrame:
    if version == "butler":
        cols = [(0, 14), (15, 28), (29, 38), (39, 45), (46, 53), (54, 62), (63, 69), (70, 74)]
        names = ["name", "t", "rv", "err", "S", "H", "counts", "exptime"]
        path = KECK.local("table1.dat.gz")
    else:
        cols = [(0, 14), (15, 28), (29, 39), (40, 47), (48, 55), (56, 63), (64, 70), (71, 75),
                (76, 82), (83, 88)]
        names = ["name", "t", "rv", "err", "S", "H", "counts", "exptime", "corr", "corr_err"]
        path = KECK_TALOR.local("table1.dat.gz")
    with gzip.open(path, "rt") as fh:
        return pd.read_fwf(fh, colspecs=cols, names=names)


def load_keck(star: Star, log: CutLog, version: str = "butler") -> pd.DataFrame:
    df = _read_keck_table(version)
    df = df[df["name"] == star.keck].copy()
    src = "Keck/HIRES" + ("" if version == "butler" else " (Tal-Or corrected)")
    log.add(src, "velocities in the table", len(df))
    df = _drop_large_errors(df, log, src)
    df["inst"] = np.where(df["t"] < HIRES_UPGRADE_BJD, "HIRES pre-2004", "HIRES post-2004")
    df["H"] = df["H"].where(df["H"] > -0.5)  # -1 marks a missing H-alpha index
    return df[["t", "rv", "err", "inst", "S", "H"]]


def load_lick(star: Star, log: CutLog, *, convert_time: bool = True) -> pd.DataFrame:
    """Lick velocities. ``convert_time=False`` keeps the table's Julian dates
    as they are (used only to test which time convention the table uses)."""
    cols = [(0, 7), (8, 19), (20, 29), (30, 36), (37, 41), (42, 44)]
    df = pd.read_fwf(LICK.local("table2.dat"), colspecs=cols,
                     names=["name", "jd", "rv", "err", "snr", "dewar"], dtype={"name": str})
    df = df[df["name"] == star.lick].copy()
    log.add("Lick/Hamilton", "velocities in the table", len(df))
    df = df.dropna(subset=["jd", "rv", "err"])
    log.add("Lick/Hamilton", "with time, velocity and error", len(df))
    jd = df["jd"].values + 2400000.0
    df["t"] = jd_utc_to_bjd_tdb(jd, star, "lick") if convert_time else jd
    df = _drop_large_errors(df, log, "Lick/Hamilton")
    early = df["t"] < LICK_EARLY_END_BJD
    df["inst"] = np.where(early, "Lick 1987-94", "Lick dewar " + df["dewar"].astype(int).astype(str))
    return df[["t", "rv", "err", "inst"]]


def load_harps(star: Star, log: CutLog, *, flare_cut: bool = False) -> pd.DataFrame:
    df = pd.read_csv(rvbank_path(star.rvbank))
    log.add("HARPS", "spectra in RVBank", len(df))
    df = df[df["f_RV"] == RVBANK_GOOD_FLAG]
    log.add("HARPS", "SERVAL flag 0 (normal spectrum)", len(df))
    df = df[(df["RV_mlc_nzp"] > -9e6) & (df["e_RV_mlc_nzp"] > 0) & np.isfinite(df["RV_mlc_nzp"])]
    log.add("HARPS", "with a velocity", len(df))
    out = pd.DataFrame({"t": df["BJD"].values, "rv": df["RV_mlc_nzp"].values,
                        "err": df["e_RV_mlc_nzp"].values})
    for col, name in HARPS_ACTIVITY.items():
        v = pd.to_numeric(df[col], errors="coerce").values.astype(float)
        out[name] = np.where(v < -9e6, np.nan, v)
    out["inst"] = np.where(out["t"] < HARPS_FIBRE_BJD, "HARPS pre-2015", "HARPS post-2015")
    if np.any(np.abs(out["t"] - HARPS_FIBRE_BJD) < 3):
        raise RuntimeError("HARPS spectra within 3 days of the fibre-exchange split; "
                           "check HARPS_FIBRE_BJD")
    out = _drop_large_errors(out, log, "HARPS")
    out = _drop_distorted_lines(out, log)
    if flare_cut:
        out = _drop_flares(out, log)
    return out


def _robust_z(v):
    med = np.nanmedian(v)
    mad = 1.4826 * np.nanmedian(np.abs(v - med))
    return (v - med) / mad if mad > 0 else np.zeros_like(v)


def _drop_distorted_lines(df: pd.DataFrame, log: CutLog) -> pd.DataFrame:
    """Remove spectra whose line shapes (SERVAL dLW, CRX) are wildly off."""
    keep = np.ones(len(df), bool)
    for inst in df["inst"].unique():
        m = (df["inst"] == inst).values
        for col in ("dLW", "CRX"):
            z = _robust_z(df.loc[m, col].values)
            keep[np.flatnonzero(m)[np.abs(z) > LINE_SHAPE_NSIGMA]] = False
    log.add("HARPS", f"line shape normal (|dLW|, |CRX| < {LINE_SHAPE_NSIGMA:.0f} robust sigma)",
            int(keep.sum()))
    return df[keep]


def _drop_flares(df: pd.DataFrame, log: CutLog) -> pd.DataFrame:
    """Remove spectra taken during flares: H-alpha far above its usual level."""
    ha = df["Halpha"]
    med = np.nanmedian(ha)
    mad = 1.4826 * np.nanmedian(np.abs(ha - med))
    keep = ~(ha > med + FLARE_HALPHA_NSIGMA * mad)
    log.add("HARPS", f"H-alpha < median + {FLARE_HALPHA_NSIGMA:.0f} sigma (no flare)", int(keep.sum()))
    return df[keep.values]


def load_elodie(star: Star, log: CutLog) -> pd.DataFrame:
    ccf = parse_elodie_ccf(elodie_ccf_path(star.elodie).read_text(encoding="latin-1"))
    log.add("ELODIE", "cross-correlation velocities in the archive", len(ccf))
    ccf = ccf[ccf["type"] == "obj"]
    log.add("ELODIE", "object fibre (not sky)", len(ccf))
    ccf = ccf[ccf["masque"] == ELODIE_MASK]
    log.add("ELODIE", f"{ELODIE_MASK} mask", len(ccf))
    ccf = ccf[ccf["imatyp"].str.startswith(ELODIE_IMATYP_PREFIX)]
    log.add("ELODIE", "simultaneous thorium reference (drift measured)", len(ccf))
    hdr = pd.read_csv(elodie_times_path(star.elodie), dtype={"datenuit": str, "imanum": str})
    df = ccf.merge(hdr, on=["datenuit", "imanum"], how="left", validate="one_to_one")
    if df["JDB1"].isna().any():
        raise RuntimeError("ELODIE spectra without header times")
    bjd_utc = df["JDB1"].astype(float) + df["JDB2"].astype(float)
    out = pd.DataFrame({
        "t": bjd_utc_to_tdb(bjd_utc.values),
        "rv": df["vfit"].values * 1000.0,  # km/s -> m/s
        "err": np.full(len(df), ELODIE_NOMINAL_ERR),
        "snr": df["sn"].values,
        "BERV": df["BERV"].astype(float).values,  # km/s, from the FITS header
        "inst": "ELODIE",
    })
    return out


def _drop_large_errors(df: pd.DataFrame, log: CutLog, source: str) -> pd.DataFrame:
    med = df["err"].median()
    keep = df["err"] <= MAX_ERROR_OVER_MEDIAN * med
    log.add(source, f"error <= {MAX_ERROR_OVER_MEDIAN:.0f}x median ({med:.1f} m/s)", int(keep.sum()))
    return df[keep]


# ---------------------------------------------------------------------------
# Assembling a star
# ---------------------------------------------------------------------------

#: Instrument display order (roughly chronological), used for colours too.
INSTRUMENT_ORDER = [
    "ELODIE", "Lick 1987-94", "Lick dewar 13", "Lick dewar 6", "Lick dewar 8",
    "HIRES pre-2004", "HIRES post-2004", "HARPS pre-2015", "HARPS post-2015",
]


def build(star_key: str, sources=("elodie", "lick", "keck", "harps"), *,
          flare_cut: bool = False, keck_version: str = "butler",
          lick_convert_time: bool = True,
          min_per_instrument: int = 5) -> tuple[RVData, CutLog]:
    """All velocities of one star from the requested sources, as an RVData.

    Instruments with fewer than ``min_per_instrument`` points are dropped:
    each costs a free zero point, so a handful of points constrain nothing.
    """
    star = STARS[star_key]
    log = CutLog()
    frames = []
    if "elodie" in sources and star.elodie:
        frames.append(load_elodie(star, log))
    if "lick" in sources and star.lick:
        frames.append(load_lick(star, log, convert_time=lick_convert_time))
    if "keck" in sources and star.keck:
        frames.append(load_keck(star, log, keck_version))
    if "harps" in sources and star.rvbank:
        frames.append(load_harps(star, log, flare_cut=flare_cut))
    df = pd.concat(frames, ignore_index=True)
    counts = df["inst"].value_counts()
    small = counts[counts < min_per_instrument].index.tolist()
    if small:
        log.add("all", f"dropped instruments with < {min_per_instrument} points: {', '.join(small)}",
                int((~df["inst"].isin(small)).sum()))
        df = df[~df["inst"].isin(small)]
    names = [n for n in INSTRUMENT_ORDER if n in set(df["inst"])]
    inst = df["inst"].map({n: i for i, n in enumerate(names)}).values
    activity = {c: df[c].values for c in df.columns
                if c not in ("t", "rv", "err", "inst") and df[c].dtype.kind == "f"}
    data = RVData(star.name, df["t"].values, df["rv"].values, df["err"].values, inst, names, activity)
    return data, log
