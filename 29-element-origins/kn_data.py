"""AT2017gfo photometry: Villar et al. (2017)'s unified table, and filters.

Villar et al. gathered 714 measurements of the kilonova from 18 papers and 46
instruments, corrected the systematic offsets they found, and flagged each
row: ``*`` used in their fit, ``X`` excluded, ``O`` a visual outlier, ``R`` a
repeated reduction of the same images, ``A`` an average of repeats. We use
exactly the rows they used (``*``, ``*,A``, ``*,R``; not ``O``), 647 of them.

Each row is tied to a transmission curve from the SVO Filter Profile Service,
as Villar et al. did. Where SVO lacks the exact instrument we use the
standard system the photometry was calibrated to (SDSS ugriz, Bessell UBVRI,
2MASS JHKs); for a blackbody that changes a magnitude by a few hundredths,
well under the 0.24 mag scatter term of the fit.
"""

from __future__ import annotations

import re

import numpy as np

from config import KN_DIR, KN_T_MERGER_MJD

PHOT_TABLE = KN_DIR / "villar2017_phot_table.dat"
FILTER_DIR = KN_DIR / "filters"

#: Instrument-specific SVO profiles: (telescope or instrument keyword, band) -> SVO ID
_SPECIFIC = {
    ("UVOT", "W2"): "Swift/UVOT.UVW2", ("UVOT", "M2"): "Swift/UVOT.UVM2",
    ("UVOT", "W1"): "Swift/UVOT.UVW1", ("UVOT", "U"): "Swift/UVOT.U",
    ("UVOT", "B"): "Swift/UVOT.B", ("UVOT", "V"): "Swift/UVOT.V",
    ("DECam", "u"): "CTIO/DECam.u", ("DECam", "g"): "CTIO/DECam.g",
    ("DECam", "r"): "CTIO/DECam.r", ("DECam", "i"): "CTIO/DECam.i",
    ("DECam", "z"): "CTIO/DECam.z", ("DECam", "Y"): "CTIO/DECam.Y",
    ("GFC", "r"): "PAN-STARRS/PS1.r", ("GFC", "i"): "PAN-STARRS/PS1.i",
    ("GFC", "z"): "PAN-STARRS/PS1.z", ("GFC", "y"): "PAN-STARRS/PS1.y",
    ("VIRCAM", "Y"): "Paranal/VISTA.Y", ("VIRCAM", "J"): "Paranal/VISTA.J",
    ("VIRCAM", "Ks"): "Paranal/VISTA.Ks",
    ("GROND", "g"): "LaSilla/GROND.g", ("GROND", "r"): "LaSilla/GROND.r",
    ("GROND", "i"): "LaSilla/GROND.i", ("GROND", "z"): "LaSilla/GROND.z",
    ("GROND", "J"): "LaSilla/GROND.J", ("GROND", "H"): "LaSilla/GROND.H",
    ("GROND", "K"): "LaSilla/GROND.K",
    ("FLAMINGOS-2", "J"): "Gemini/Flamingos2.J", ("FLAMINGOS-2", "H"): "Gemini/Flamingos2.H",
    ("FLAMINGOS-2", "Ks"): "Gemini/Flamingos2.Ks",
    ("GMOS", "g"): "Gemini/GMOS-S.g", ("GMOS", "r"): "Gemini/GMOS-S.r",
    ("GMOS", "i"): "Gemini/GMOS-S.i", ("GMOS", "z"): "Gemini/GMOS-S.z",
    ("HAWKI", "Ks"): "Paranal/HAWKI.Ks",
    ("FourStar", "J1"): "LCO/FourStar.J1",
    ("VISIR", "J8.9"): "Paranal/VISIR.J89",
    ("HSC", "z"): "Subaru/HSC.z",
    ("Skymapper", "g"): "SkyMapper/SkyMapper.g", ("Skymapper", "r"): "SkyMapper/SkyMapper.r",
    ("Skymapper", "i"): "SkyMapper/SkyMapper.i",
}
#: Standard systems for everything else.
_GENERIC = {
    "u": "SLOAN/SDSS.u", "g": "SLOAN/SDSS.g", "r": "SLOAN/SDSS.r", "i": "SLOAN/SDSS.i",
    "z": "SLOAN/SDSS.z", "y": "PAN-STARRS/PS1.y", "Y": "Paranal/VISTA.Y",
    "U": "Generic/Bessell.U", "B": "Generic/Bessell.B", "V": "Generic/Bessell.V",
    "R": "Generic/Bessell.R", "I": "Generic/Bessell.I",
    "J": "2MASS/2MASS.J", "H": "2MASS/2MASS.H", "Ks": "2MASS/2MASS.Ks", "K": "2MASS/2MASS.Ks",
    # MASTER's unfiltered "W" band is calibrated as 0.2B + 0.8R (Lipunov et al.
    # 2017); R is the closest standard band.
    "W": "Generic/Bessell.R",
}


def svo_id(telescope: str, instrument: str, band: str) -> str:
    if telescope == "HST":
        band = "F850LP" if band == "F850W" else band
        if instrument.startswith("ACS"):
            return f"HST/ACS_WFC.{band}"
        if instrument == "WFC3/IR":
            return f"HST/WFC3_IR.{band}"
        chip = "WFC3_UVIS1" if instrument.endswith("UVIS1") else "WFC3_UVIS2"
        return f"HST/{chip}.{band}"
    for key in (instrument, telescope):
        if (key, band) in _SPECIFIC:
            return _SPECIFIC[(key, band)]
    return _GENERIC.get(band, "")


def _clean(cell: str) -> str:
    return re.sub(r"\s+", " ", cell.replace("\\\\", "").replace("$", "")).strip()


def read_villar_table(path=PHOT_TABLE) -> dict:
    """Parse Table 3 of Villar et al. (2017) from its LaTeX source."""
    rows = []
    for line in path.read_text().splitlines():
        if not re.match(r"\s*\d{5}\.\d+\s*&", line):
            continue
        # "B\&C" (a telescope name) contains an escaped ampersand.
        cells = [_clean(c) for c in line.replace("\\&", "and").split("&")]
        if len(cells) < 10:
            raise ValueError(f"unparsed row: {line!r}")
        mjd, phase, inst, tel, band, mag, err, dmag, ref, note = cells[:10]
        upper = mag.startswith(">")
        rows.append(dict(
            mjd=float(mjd), phase=float(phase), instrument=inst, telescope=tel, band=band,
            mag=float(mag.lstrip(">")), err=float(err) if re.match(r"^[\d.]+$", err) else np.nan,
            upper=upper, ref=ref.replace("citeauthor", "").strip("{}\\ "),
            note=note.replace(" ", "").upper(),
        ))
    out = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    flags = [set(n.split(",")) for n in out["note"]]
    out["used"] = np.array(["*" in f and "O" not in f and "X" not in f for f in flags])
    out["t"] = out["mjd"] - KN_T_MERGER_MJD
    out["svo"] = np.array([svo_id(t, i, b) for t, i, b in
                           zip(out["telescope"], out["instrument"], out["band"])])
    unmapped = out["used"] & (out["svo"] == "")
    if unmapped.any():
        raise ValueError(f"no filter for used rows: {set(out['band'][unmapped])}")
    return out


def filter_path(svo: str):
    return FILTER_DIR / (svo.replace("/", "__") + ".dat")


def read_filter(svo: str):
    """Wavelength (Angstrom) and transmission of an SVO filter."""
    lam, tr = np.loadtxt(filter_path(svo), unpack=True)
    keep = tr > 0
    # Keep one zero on each side so the edges integrate cleanly.
    idx = np.flatnonzero(keep)
    lo, hi = max(idx[0] - 1, 0), min(idx[-1] + 2, lam.size)
    return lam[lo:hi], np.clip(tr[lo:hi], 0.0, None)
