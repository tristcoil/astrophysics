"""Loaders for every light curve in the project, as plain numpy arrays.

Times are kept in each source's own convention and converted here to one:
NGC 5548 epochs are THJD = HJD - 2,450,000 (AGN STORM's convention); the AGN
Watch years are JD - 2,400,000 (their convention, kept because the published
season boundaries use it); SDSS-RM uses MJD.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

import config as C

VIZ = C.DATA_DIR / "vizier"


@dataclass
class LightCurve:
    t: np.ndarray
    y: np.ndarray
    e: np.ndarray
    label: str = ""
    unit: str = ""

    def __post_init__(self):
        o = np.argsort(self.t, kind="stable")
        self.t, self.y, self.e = (np.asarray(a, float)[o] for a in (self.t, self.y, self.e))

    def __len__(self):
        return len(self.t)

    def window(self, lo=-np.inf, hi=np.inf) -> "LightCurve":
        m = (self.t >= lo) & (self.t <= hi)
        return LightCurve(self.t[m], self.y[m], self.e[m], self.label, self.unit)

    def shifted(self, dt: float) -> "LightCurve":
        return LightCurve(self.t + dt, self.y, self.e, self.label, self.unit)


# ---------------------------------------------------------------------------
# AGN STORM 2014
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _pei_rows():
    rows = []
    for line in (VIZ / "J_ApJ_837_131" / "table2.dat").read_text().splitlines():
        if not line.strip():
            continue

        def num(a, b):
            s = line[a - 1:b].strip()
            return float(s) if s else np.nan

        rows.append((num(1, 7), line[8:14].strip(), num(16, 21), num(23, 27), num(29, 35), num(37, 41),
                     num(43, 49), num(51, 55), num(57, 63), num(65, 69)))
    return rows


def storm_spectroscopy(telescope: str | None = "MDM", truncate: bool = True) -> dict:
    """Pei et al. (2017) ground-based spectroscopy of NGC 5548 (Table 2).

    Returns light curves for the 5100 A continuum (1e-15 erg/s/cm^2/A), H-beta
    from a linear continuum ("Hb"), H-beta and He II 4686 from spectral
    decomposition ("HbSD", "HeIISD"), all in 1e-15 erg/s/cm^2 and including the
    narrow lines. ``telescope=None`` keeps all five sites (already
    intercalibrated by Pei et al. for 5100 A and Hb; the SD fluxes exist for MDM
    only). ``truncate`` drops the last 10 MDM epochs (THJD > 6828.75) with
    inconsistent calibration, as Pei et al. did for their lags.
    """
    rows = [r for r in _pei_rows() if telescope is None or r[1] == telescope]
    if truncate:
        rows = [r for r in rows if r[0] <= C.STORM_TRUNCATE_THJD]
    a = np.array([[r[0], *r[2:]] for r in rows], float)
    t = a[:, 0]
    out = {
        "F5100": LightCurve(t, a[:, 1], a[:, 2], "5100 Å continuum", "1e-15 erg/s/cm²/Å"),
        "Hb": LightCurve(t, a[:, 3], a[:, 4], "Hβ", "1e-15 erg/s/cm²"),
    }
    sd = np.isfinite(a[:, 5])
    out["HbSD"] = LightCurve(t[sd], a[sd, 5], a[sd, 6], "Hβ (spectral decomposition)", "1e-15 erg/s/cm²")
    out["HeIISD"] = LightCurve(t[sd], a[sd, 7], a[sd, 8], "He II 4686", "1e-15 erg/s/cm²")
    return out


@lru_cache(maxsize=None)
def _fausnaugh_hst():
    a = np.loadtxt(VIZ / "J_ApJ_821_56" / "table4.dat")
    return a


def storm_hst_continuum(wavelength: float = 1157.5) -> LightCurve:
    """HST/COS continuum flux density (Fausnaugh et al. 2016, Table 4).

    ``wavelength`` in {1157.5, 1367, 1478.5, 1746} A; flux in 1e-15 erg/s/cm^2/A.
    """
    a = _fausnaugh_hst()
    m = np.isclose(a[:, 0], wavelength)
    if not m.any():
        raise ValueError(f"no HST continuum at {wavelength} A")
    return LightCurve(a[m, 1] - 50000.0, a[m, 2], a[m, 3], f"{wavelength:g} Å continuum (HST)",
                      "1e-15 erg/s/cm²/Å")


def storm_photometry() -> dict[str, LightCurve]:
    """Ground-based broad-band continuum (Fausnaugh et al. 2016, Table 3).

    One light curve per filter (u g r i z B V R I), all telescopes merged as
    published (intercalibrated by the authors); flux in 1e-15 erg/s/cm^2/A.
    """
    out: dict[str, list] = {}
    for line in (VIZ / "J_ApJ_821_56" / "table3.dat").read_text().splitlines():
        if not line.strip():
            continue
        f = line[0]
        out.setdefault(f, []).append((float(line[2:13]) - 50000.0, float(line[14:21]), float(line[22:28])))
    res = {}
    for f, rows in out.items():
        a = np.array(rows)
        res[f] = LightCurve(a[:, 0], a[:, 1], a[:, 2], f"{f} band", "1e-15 erg/s/cm²/Å")
    return res


def storm_uv_lines() -> dict[str, LightCurve]:
    """HST/COS broad UV lines and the 1367 A continuum (De Rosa et al. 2015).

    Line fluxes in 1e-13 erg/s/cm^2 (narrow lines included), the continuum in
    1e-15 erg/s/cm^2/A. C IV and He II come from the G160M grating, observed
    ~1 hour after the G130M exposure that gives the others.
    """
    a = np.loadtxt(VIZ / "J_ApJ_806_128" / "table2.dat")
    t130, t160 = a[:, 0], a[:, 9]
    return {
        "F1367": LightCurve(t130, a[:, 1], a[:, 2], "1367 Å continuum", "1e-15 erg/s/cm²/Å"),
        "Lya": LightCurve(t130, a[:, 3], a[:, 4], "Lyα", "1e-13 erg/s/cm²"),
        "NV": LightCurve(t130, a[:, 5], a[:, 6], "N V", "1e-13 erg/s/cm²"),
        "SiIV": LightCurve(t130, a[:, 7], a[:, 8], "Si IV + O IV]", "1e-13 erg/s/cm²"),
        "CIV": LightCurve(t160, a[:, 10], a[:, 11], "C IV", "1e-13 erg/s/cm²"),
        "HeII": LightCurve(t160, a[:, 12], a[:, 13], "He II 1640 + O III]", "1e-13 erg/s/cm²"),
    }


# ---------------------------------------------------------------------------
# AGN Watch and later campaigns, 1988-2008 (Peterson et al. 2013)
# ---------------------------------------------------------------------------


def watch() -> dict[str, LightCurve]:
    """NGC 5548's host-subtracted 5100 A AGN continuum and broad H-beta, 1988-2008.

    JD - 2,400,000. Continuum in 1e-15 erg/s/cm^2/A (observed frame, starlight
    removed), H-beta in 1e-13 erg/s/cm^2 (narrow component removed).
    """
    c = np.loadtxt(VIZ / "J_ApJ_779_109" / "table5.dat")
    h = np.loadtxt(VIZ / "J_ApJ_779_109" / "table6.dat")
    return {
        "F5100": LightCurve(c[:, 0], c[:, 1], c[:, 2], "5100 Å AGN continuum", "1e-15 erg/s/cm²/Å"),
        "Hb": LightCurve(h[:, 0], h[:, 1], h[:, 2], "broad Hβ", "1e-13 erg/s/cm²"),
    }


# ---------------------------------------------------------------------------
# SDSS-RM
# ---------------------------------------------------------------------------


def sdssrm(rmid: int, line: str) -> LightCurve | None:
    """One SDSS-RM light curve: ``line`` in {cont, hb, ha, mg2, c4}; MJD.

    The continuum is the PyCALI-merged 11-year light curve (f_nu, 1e-28
    erg/s/cm^2/Hz) from SDSS, BOK, CFHT and Pan-STARRS photometry; lines are
    PrepSpec fluxes (1e-17 erg/s/cm^2) from the 2014-2020 spectroscopy.
    """
    p = C.DATA_DIR / "sdssrm" / f"rm{rmid:03d}" / f"{line}_lc.txt"
    if not p.exists():
        return None
    a = np.loadtxt(p, usecols=(0, 1, 2), comments="#", ndmin=2)
    good = np.isfinite(a).all(axis=1) & (a[:, 2] > 0)
    a = a[good]
    if len(a) == 0:
        return None
    # A handful of continuum epochs repeat the same MJD from different facilities:
    # PyCCF needs strictly increasing times, so equal times are averaged.
    t, inv = np.unique(a[:, 0], return_inverse=True)
    if len(t) < len(a):
        w = 1.0 / a[:, 2] ** 2
        sw = np.bincount(inv, w)
        y = np.bincount(inv, w * a[:, 1]) / sw
        e = 1.0 / np.sqrt(sw)
        return LightCurve(t, y, e, line, "")
    return LightCurve(a[:, 0], a[:, 1], a[:, 2], line, "")
