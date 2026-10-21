"""Filters, magnitudes and fluxes; loaders that put every photometry source into one table.

Each loader returns a pandas DataFrame with the columns

    mjd, phase, band, mag, err, system, source

``band`` is a canonical filter name (see ``BANDS``), ``system`` is "Vega" or
"AB", ``phase`` is days since the explosion epoch in ``config.SUPERNOVAE``.
"""

from __future__ import annotations

import io
import json
import re
import tarfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

import config as C

C_ANG = C.C_CGS * 1e8  # speed of light in Angstrom/s
AB_ZP_JY = 3631.0

#: Canonical band name -> SVO filter ID. Johnson-Cousins bands use Bessell's
#: (1990) realisation; JHK photometry calibrated to 2MASS uses 2MASS curves.
BANDS = {
    "UVW2": "Swift/UVOT.UVW2", "UVM2": "Swift/UVOT.UVM2", "UVW1": "Swift/UVOT.UVW1",
    "U": "Generic/Bessell.U", "B": "Generic/Bessell.B", "V": "Generic/Bessell.V",
    "R": "Generic/Bessell.R", "I": "Generic/Bessell.I",
    "u": "SLOAN/SDSS.u", "g": "SLOAN/SDSS.g", "r": "SLOAN/SDSS.r", "i": "SLOAN/SDSS.i",
    "z": "SLOAN/SDSS.z",
    "J": "2MASS/2MASS.J", "H": "2MASS/2MASS.H", "Ks": "2MASS/2MASS.Ks",
    "ztf_g": "Palomar/ZTF.g", "ztf_r": "Palomar/ZTF.r",
}


@dataclass(frozen=True)
class Filter:
    name: str
    svo_id: str
    pivot: float  # Angstrom
    zp_vega_jy: float  # flux density of a Vega-system zero magnitude, Jy
    wave: np.ndarray  # Angstrom
    trans: np.ndarray


def _parse_svo(path) -> tuple[dict, np.ndarray, np.ndarray]:
    root = ET.parse(path).getroot()
    params = {}
    for el in root.iter():
        if el.tag.endswith("PARAM") and "name" in el.attrib:
            params[el.attrib["name"]] = el.attrib.get("value")
    rows = []
    for tr in root.iter():
        if tr.tag.endswith("TR"):
            vals = [td.text for td in tr if td.tag.endswith("TD")]
            if len(vals) >= 2:
                rows.append((float(vals[0]), float(vals[1])))
    arr = np.array(rows)
    return params, arr[:, 0], arr[:, 1]


@lru_cache(maxsize=None)
def get_filter(name: str) -> Filter:
    svo_id = BANDS[name]
    path = C.DATA_DIR / "filters" / (svo_id.replace("/", "_") + ".xml")
    params, wave, trans = _parse_svo(path)
    if params.get("MagSys") != "Vega":
        raise ValueError(f"{svo_id}: expected Vega zero point, got {params.get('MagSys')}")
    return Filter(name, svo_id, float(params["WavelengthPivot"]), float(params["ZeroPoint"]),
                  wave, trans)


def pivot(band: str) -> float:
    return get_filter(band).pivot


def mag_to_flam(mag, band: str, system: str):
    """Magnitude -> flux density f_lambda (erg/s/cm^2/A) at the band's pivot wavelength."""
    f = get_filter(band)
    zp = AB_ZP_JY if system == "AB" else f.zp_vega_jy
    fnu = zp * 1e-23 * 10 ** (-0.4 * np.asarray(mag, float))  # erg/s/cm^2/Hz
    return fnu * C_ANG / f.pivot**2


def synthetic_mag(wave, flam, band: str, system: str = "AB") -> float:
    """Synthetic magnitude of a spectrum through a filter (photon-counting).

    m = -2.5 log10( Int f_lambda T lambda dlambda / Int f_nu,ref T c/lambda dlambda ),
    with f_nu,ref the AB (or Vega-system) zero point.
    """
    f = get_filter(band)
    t = np.interp(wave, f.wave, f.trans, left=0.0, right=0.0)
    if np.any(t[[0, -1]] > 0.01 * t.max()):
        raise ValueError(f"spectrum does not cover the {band} filter")
    num = np.trapezoid(flam * t * wave, wave)
    zp = AB_ZP_JY if system == "AB" else f.zp_vega_jy
    den = np.trapezoid(zp * 1e-23 * C_ANG / wave**2 * t * wave, wave)
    return float(-2.5 * np.log10(num / den))


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

COLUMNS = ["mjd", "phase", "band", "mag", "err", "system", "source"]

#: Errors assigned where a source gives none: photoelectric UBVRI of SN 1987A
#: from SAAO and CTIO is good to ~0.02 mag; other unflagged entries 0.05 mag.
DEFAULT_ERR = {"SN1987A": 0.02}


def _finish(df: pd.DataFrame, sn: str) -> pd.DataFrame:
    df = df.copy()
    df["phase"] = df["mjd"] - C.SUPERNOVAE[sn]["t0_mjd"]
    df = df.dropna(subset=["mag"])
    # Drop exact duplicates (the same point reached the catalog via two papers).
    key = pd.DataFrame({"t": df["mjd"].round(2), "b": df["band"], "m": df["mag"].round(2)})
    df = df[~key.duplicated()]
    return df[COLUMNS].sort_values(["band", "mjd"]).reset_index(drop=True)


def load_osc(sn: str) -> pd.DataFrame:
    """Photometry from the Open Supernova Catalog JSON.

    Keeps Johnson-Cousins UBVRI (ground) and the Swift UVOT ultraviolet filters.
    UVOT's own u, b, v filters are not Johnson's and are skipped, as are upper
    limits, unfiltered data, host-included or model points.
    """
    d = json.loads((C.DATA_DIR / "osc" / f"{sn}.json").read_text())[sn]
    names = {s["alias"]: s.get("reference", s.get("name", "")) for s in d["sources"]}
    rows = []
    for p in d["photometry"]:
        if "magnitude" not in p or p.get("upperlimit") or p.get("includeshost") or p.get("model"):
            continue
        band = p.get("band")
        swift = p.get("telescope") == "Swift" or p.get("instrument") == "UVOT"
        if swift:
            if band not in ("UVW2", "UVM2", "UVW1") or p.get("system") != "Vega":
                continue  # the copies without a system tag duplicate the Vega ones
        else:
            band = {"Rc": "R", "Ic": "I"}.get(band, band)
            if band not in ("U", "B", "V", "R", "I"):
                continue
        t = p["time"][0] if isinstance(p["time"], list) else p["time"]
        err = p.get("e_magnitude") or p.get("e_upper_magnitude")
        err = float(err) if err not in (None, "") else DEFAULT_ERR.get(sn, 0.05)
        src = ";".join(names.get(a, a) for a in p.get("source", "").split(","))
        rows.append((float(t), band, float(p["magnitude"]), err, "Vega", src))
    df = pd.DataFrame(rows, columns=["mjd", "band", "mag", "err", "system", "source"])
    return _finish(df, sn)


def _arxiv_member(label_prefix: str, suffix: str) -> str:
    path = next((C.DATA_DIR / "arxiv").glob(f"{label_prefix}_*.tar.gz"))
    with tarfile.open(path) as tar:
        for m in tar.getmembers():
            if m.name.endswith(suffix):
                return tar.extractfile(m).read().decode("latin-1")
    raise FileNotFoundError(f"{suffix} not in {path.name}")


def load_matheson2012() -> pd.DataFrame:
    """SN 2011fe JHKs (Matheson et al. 2012, Table 1), on the 2MASS (Vega) system."""
    tex = _arxiv_member("matheson2012", "astroph.tex")
    block = tex.split("caption{Apparent magnitudes")[1].split(r"\enddata")[0]
    rows = []
    for line in block.splitlines():
        m = re.match(r"\s*(5\d{4}\.\d+)\s*&(.*)\\\\", line)
        if not m:
            continue
        cells = [c.strip() for c in m.group(2).split("&")]
        for band, cell in zip(("J", "H", "Ks"), cells):
            mm = re.match(r"([\d.]+)\s*\((\d+)\)", cell)
            if mm:
                rows.append((float(m.group(1)), band, float(mm.group(1)), int(mm.group(2)) / 1000.0,
                             "Vega", "Matheson et al. (2012)"))
    return _finish(pd.DataFrame(rows, columns=["mjd", "band", "mag", "err", "system", "source"]),
                   "SN2011fe")


def load_amanullah2014() -> pd.DataFrame:
    """SN 2014J JHK (Mount Abu) from Amanullah et al. (2014), natural Vega system."""
    path = C.DATA_DIR / "vizier" / "J_ApJ_788_L21" / "table1.dat"
    rows = []
    for line in path.read_text().splitlines():
        if len(line) < 29:
            continue
        band = line[13:18].strip()
        if band in ("J", "H", "K"):
            rows.append((float(line[0:7]), {"K": "Ks"}[band] if band == "K" else band,
                         float(line[19:24]), float(line[25:29]), "Vega", "Amanullah et al. (2014)"))
    return _finish(pd.DataFrame(rows, columns=["mjd", "band", "mag", "err", "system", "source"]),
                   "SN2014J")


def load_shappee2017() -> pd.DataFrame:
    """Late SN 2011fe photometry from Shappee et al. (2017): LBT B, V, R and J."""
    path = C.DATA_DIR / "vizier" / "J_ApJ_841_48" / "table1.dat"
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        band = parts[1]
        if band in ("B", "V", "R", "J"):
            rows.append((float(parts[0]) - 2400000.5, band, float(parts[2]), float(parts[3]),
                         "Vega", "Shappee et al. (2017)"))
    return _finish(pd.DataFrame(rows, columns=["mjd", "band", "mag", "err", "system", "source"]),
                   "SN2011fe")


def load_hsu2025() -> pd.DataFrame:
    """SN 2023ixf, Las Cumbres Observatory UBVgriz (Hsu et al. 2025), AB magnitudes."""
    rows = []
    for line in (C.DATA_DIR / "sn2023ixf" / "hsu2025_dbf1.txt").read_text().splitlines():
        if not line[:1].isdigit():
            continue
        if line[41:46].strip() == "True":  # non-detection flag
            continue
        rows.append((float(line[0:12]), line[13:14], float(line[15:21]), float(line[22:27]),
                     "AB", "Hsu et al. (2025)"))
    return _finish(pd.DataFrame(rows, columns=["mjd", "band", "mag", "err", "system", "source"]),
                   "SN2023ixf")


def load_singh2024() -> pd.DataFrame:
    """SN 2023ixf from Singh et al. (2024): Swift UV, optical and JHKs, AB magnitudes.

    Swift's own u, b, v filters are skipped (ground UBV covers those wavelengths).
    """
    rows = []
    for line in (C.DATA_DIR / "sn2023ixf" / "singh2024_dbf2.txt").read_text().splitlines():
        if not line.startswith("246"):
            continue
        band, tel = line[29:33].strip(), line[34:47].strip()
        if tel == "SWIFT-UVOT" and band not in ("UVW2", "UVM2", "UVW1"):
            continue
        err = max(float(line[24:28]), 0.01)  # a few entries quote 0.00
        rows.append((float(line[0:10]) - 2400000.5, band, float(line[18:23]), err, "AB",
                     f"Singh et al. (2024) {tel}"))
    return _finish(pd.DataFrame(rows, columns=["mjd", "band", "mag", "err", "system", "source"]),
                   "SN2023ixf")


def load_alerce(oid: str, sn: str | None = None) -> pd.DataFrame:
    """ZTF alert photometry (difference-image PSF magnitudes, AB) from ALeRCE."""
    d = json.loads((C.DATA_DIR / "ztf" / "lightcurves" / f"{oid}.json").read_text())
    rows = []
    for p in d["detections"]:
        if p.get("isdiffpos") not in (1, "1", "t", True) or p.get("magpsf") is None:
            continue
        band = {1: "ztf_g", 2: "ztf_r"}.get(p["fid"])
        if band:
            rows.append((p["mjd"], band, p["magpsf"], p["sigmapsf"], "AB", "ZTF/ALeRCE"))
    df = pd.DataFrame(rows, columns=["mjd", "band", "mag", "err", "system", "source"])
    if sn is None:
        df["phase"] = np.nan
        return df[COLUMNS].sort_values(["band", "mjd"]).reset_index(drop=True)
    return _finish(df, sn)


#: Data-quality cuts, each justified by comparison with independent data.
#: SN 2011fe: the KAIT/Nickel photometry (Silverman et al. 2012, re-released by
#: Stahl et al. 2019) agrees with the other sources to ~0.05 mag until day
#: ~150, then drifts faint by 0.3 mag (day 200-250) to 1.2-1.4 mag (day
#: 300-360). The other series (Munari et al. 2013; Tsvetkov et al. 2013) join
#: smoothly onto the independent LBT photometry of Shappee et al. (2017) at day
#: 506-653, so the drift is in the KAIT/Nickel late-time data.
LATE_CUTS = {"SN2011fe": [("Stahl et al. (2019)", 150.0)]}


def load_supernova(sn: str) -> pd.DataFrame:
    """All photometry used to build the bolometric light curve of one supernova."""
    if sn == "SN2011fe":
        osc = load_osc(sn)
        for src, t_max in LATE_CUTS[sn]:
            osc = osc[~(osc["source"].str.contains(src, regex=False) & (osc["phase"] > t_max))]
        parts = [osc, load_matheson2012(), load_shappee2017()]
    elif sn == "SN2014J":
        parts = [load_osc(sn), load_amanullah2014()]
    elif sn == "SN1987A":
        parts = [load_osc(sn)]
    elif sn == "SN2023ixf":
        parts = [load_hsu2025(), load_singh2024()]
    else:
        raise KeyError(sn)
    return pd.concat(parts, ignore_index=True).sort_values(["band", "mjd"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Published bolometric light curves of SN 1987A
# ---------------------------------------------------------------------------


def load_sn1987a_uvoir() -> pd.DataFrame:
    """UVOIR luminosity of SN 1987A, days 14-935 (Suntzeff et al. 1991; Bouchet et al. 1991).

    As tabulated by D. Jeffery; days since the neutrino burst, log10 L (erg/s)
    for 50 kpc. ``table`` 5-6 are Suntzeff et al., 3-4 Bouchet et al.
    """
    text = (C.DATA_DIR / "sn1987a" / "jeffery_lc_bol2.dat").read_text()
    rows = []
    for line in text.split("END:")[1].splitlines():
        parts = line.split()
        if len(parts) >= 2:
            table = int(parts[2]) if len(parts) > 2 else 5
            src = "Suntzeff et al. (1991)" if table in (5, 6) else "Bouchet et al. (1991)"
            rows.append((float(parts[0]), float(parts[1]), table, src))
    return pd.DataFrame(rows, columns=["phase", "log_l", "table", "source"])


def load_sn1987a_late_uvoir() -> pd.DataFrame:
    """Late UVOIR points of SN 1987A (days 282-1854) read from Seitenzahl et al. 2014, Fig. 2.

    The figure is a matplotlib EPS: the 28 red hexagons are the UVOIR data of
    Suntzeff et al. (1992, 1997). Positions are converted to data units with the
    axis tick marks (0-4000 days; log L = 34-42), which recovers the published
    days as integers, i.e. the extraction is exact to the plotting precision.
    """
    path = next((C.DATA_DIR / "arxiv").glob("seitenzahl2014_*.tar.gz"))
    with tarfile.open(path) as tar:
        eps = tar.extractfile("f2.eps").read().decode("latin-1").splitlines()
    main_clip = (446.4, 345.6, 72.0, 43.2)
    clip, pts, xticks, yticks = None, [], [], []
    for line in eps:
        s = line.strip()
        m = re.match(r"^([\d.\-e]+) ([\d.\-e]+) ([\d.\-e]+) ([\d.\-e]+) clipbox$", s)
        if m:
            clip = tuple(float(v) for v in m.groups())
            continue
        m = re.match(r"^([\d.\-e]+) ([\d.\-e]+) p0_0$", s)  # red hexagon marker
        if m and clip == main_clip:
            pts.append((float(m.group(1)), float(m.group(2))))
            continue
        m = re.match(r"^([\d.\-e]+) ([\d.\-e]+) o$", s)  # tick marks
        if m:
            x, y = float(m.group(1)), float(m.group(2))
            if y == 43.2 and 72.0 < x < 518.4:  # (the corners are y-axis ticks)
                xticks.append(x)
            if x == 72.0 and 43.2 <= y <= 388.8:
                yticks.append(y)
    xt = np.unique(np.round(xticks, 3))  # ticks at 0, 1000, ..., 4000 days
    yt = np.unique(np.round(yticks, 3))  # ticks at log L = 34, 35, ..., 42
    if len(xt) != 5 or len(yt) != 9:
        raise ValueError(f"unexpected tick layout in Seitenzahl et al. Fig. 2: {xt}, {yt}")
    px = np.array(pts)
    days = (px[:, 0] - xt[0]) / (xt[-1] - xt[0]) * 4000.0
    log_l = 34.0 + (px[:, 1] - yt[0]) / (yt[-1] - yt[0]) * 8.0
    df = pd.DataFrame({"phase": np.round(days, 1), "log_l": np.round(log_l, 3),
                       "source": "Suntzeff et al. (1992, 1997) via Seitenzahl et al. (2014)"})
    return df.sort_values("phase").reset_index(drop=True)
