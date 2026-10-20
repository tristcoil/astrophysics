"""Parse the AGN Black Hole Mass Database (Bentz & Katz 2015) into tables.

Each AGN has a page listing every published reverberation campaign: the line,
the Julian-date range monitored, the lags (JAVELIN, ICCF centroid and peak,
all in the host's rest frame), the rms line width (sigma_line and FWHM), and
the host-corrected 5100 A AGN luminosity where one was measured. The database
does not store stellar velocity dispersions.
"""

from __future__ import annotations

import re
from functools import lru_cache

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

import config as C

FOLDER = C.DATA_DIR / "agnmass"

#: Wavelength tag on the database pages -> line name.
LINES = {
    "1215": "Lya", "1240": "NV", "1400": "SiIV", "1549": "CIV", "1640": "HeII1640", "1909": "CIII]",
    "2798": "MgII", "4340": "Hg", "4686": "HeII", "4861": "Hb", "5876": "HeI", "6563": "Ha",
    "4471": "HeI4471", "1661": "OIII]", "3727": "[OII]", "4102": "Hd",
}

_NUM = r"[-+]?\d+(?:\.\d+)?"


def _asym(text: str):
    """'5.8 ( + 0.4 / - 0.4 )' -> (5.8, 0.4 minus, 0.4 plus); '...' -> NaNs."""
    m = re.match(rf"\s*({_NUM})\s*\(\s*\+\s*({_NUM})\s*/\s*-\s*({_NUM})\s*\)", text)
    if not m:
        m1 = re.match(rf"\s*({_NUM})\s*$", text)
        if m1:
            return float(m1.group(1)), np.nan, np.nan
        return np.nan, np.nan, np.nan
    v, plus, minus = (float(x) for x in m.groups())
    return v, minus, plus


def _sym(text: str):
    """'2576 (+/- 389)' -> (2576, 389)."""
    m = re.match(rf"\s*({_NUM})\s*\(\s*\+/-\s*({_NUM})\s*\)", text)
    if not m:
        return np.nan, np.nan
    return float(m.group(1)), float(m.group(2))


def _page(path):
    soup = BeautifulSoup(path.read_text(encoding="latin-1"), "lxml")
    text = soup.get_text(" ", strip=True)
    name = text.split("Alternate Names")[0].strip().split()[-1]
    z = float(re.search(rf"z\s*=\s*({_NUM})", text).group(1))
    dl = re.search(rf"D\s*L\s*=\s*({_NUM})\s*Mpc", text)
    meas = []
    body = soup.find("table", class_="body")
    if body is None:
        return name, z, (float(dl.group(1)) if dl else np.nan), meas
    for tr in body.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) != 11:
            continue
        cells = [td.get_text(" ", strip=True) for td in tds]
        wl = re.search(r"λ\s*(\d+)", cells[0])
        line = LINES.get(wl.group(1), "λ" + wl.group(1)) if wl else cells[0]
        jd = re.findall(r"\d+(?:\.\d+)?", cells[1])
        jav = _asym(cells[2])
        cen = _asym(cells[4])
        pk = _asym(cells[5])
        sig = _sym(cells[6])
        fw = _sym(cells[7])
        logl = _sym(cells[9])
        meas.append(dict(
            object=name, line=line, jd_start=float(jd[0]) if jd else np.nan,
            jd_end=float(jd[1]) if len(jd) > 1 else np.nan,
            tau_jav=jav[0], tau_jav_lo=jav[1], tau_jav_hi=jav[2], ref_jav=cells[3],
            tau_cent=cen[0], tau_cent_lo=cen[1], tau_cent_hi=cen[2],
            tau_peak=pk[0], tau_peak_lo=pk[1], tau_peak_hi=pk[2],
            sigma_line=sig[0], sigma_line_err=sig[1], fwhm=fw[0], fwhm_err=fw[1],
            ref_rm=cells[8], logL5100=logl[0], logL5100_err=logl[1], ref_L=cells[10],
        ))
    return name, z, (float(dl.group(1)) if dl else np.nan), meas


@lru_cache(maxsize=None)
def measurements() -> pd.DataFrame:
    """Every reverberation measurement in the database, one row per line per campaign."""
    rows, objs = [], []
    for p in sorted(FOLDER.glob("details_*.html")):
        name, z, dl, meas = _page(p)
        objs.append(dict(object=name, z=z, DL_Mpc=dl, page=p.name))
        for m in meas:
            m["z"] = z
            rows.append(m)
    df = pd.DataFrame(rows)
    df["year"] = df["jd_start"].map(_jd_year)
    return df


@lru_cache(maxsize=None)
def objects() -> pd.DataFrame:
    out = []
    for p in sorted(FOLDER.glob("details_*.html")):
        name, z, dl, meas = _page(p)
        out.append(dict(object=name, z=z, DL_Mpc=dl, n_meas=len(meas)))
    return pd.DataFrame(out)


def _jd_year(jd):
    """Calendar year of a (JD - 2,400,000) date, or of an MJD; NaN-safe."""
    if not np.isfinite(jd):
        return np.nan
    return 1858.87 + jd / 365.25


def hbeta_rl() -> pd.DataFrame:
    """H-beta campaigns with both a centroid lag and a host-corrected L5100."""
    d = measurements()
    d = d[(d.line == "Hb") & np.isfinite(d.tau_cent) & np.isfinite(d.logL5100) & (d.tau_cent > 0)].copy()
    return d.reset_index(drop=True)
