"""Readers for the OGLE Collection of Variable Stars (OCVS) files.

The catalogue files are fixed-width text (formats in each directory's README);
missing values are written as ``-``. Multi-mode Cepheids list their periods
from the longest (lowest-order mode) to the shortest.
"""

from __future__ import annotations

import re
from functools import lru_cache

import numpy as np

import config as C

#: Pulsation-mode files and the modes of their period columns, longest first.
MODE_FILES = {
    "cepF.dat": ("F",),
    "cep1O.dat": ("1O",),
    "cep2O.dat": ("2O",),
    "cepF1O.dat": ("F", "1O"),
    "cep1O2O.dat": ("1O", "2O"),
    "cep1O3O.dat": ("1O", "3O"),
    "cep2O3O.dat": ("2O", "3O"),
    "cepF1O2O.dat": ("F", "1O", "2O"),
    "cep1O2O3O.dat": ("1O", "2O", "3O"),
}

#: Byte ranges (0-based, end-exclusive) of one period block in the mode files.
#: The first block starts at byte 34 (1-based), later ones every 66 bytes.
_BLOCK = {"P": (33, 44), "eP": (45, 54), "T0": (56, 66), "ampI": (68, 73),
          "R21": (75, 80), "phi21": (81, 86), "R31": (88, 93), "phi31": (94, 99)}
_BLOCK2 = {"P": (100, 110), "eP": (111, 120), "T0": (122, 132), "ampI": (134, 139),
           "R21": (141, 146), "phi21": (147, 152), "R31": (154, 159), "phi31": (160, 165)}
_BLOCK3 = {k: (a + 66, b + 66) for k, (a, b) in _BLOCK2.items()}


def _num(s: str) -> float:
    s = s.strip()
    return float(s) if s and s != "-" else np.nan


def _dms(sign: str, d: str, m: str, s: str) -> float:
    v = int(d) + int(m) / 60 + float(s) / 3600
    return -v if sign == "-" else v


@lru_cache(maxsize=None)
def read_ident(gal: str) -> dict:
    """ident.dat: OGLE ID, mode string, J2000 position, survey IDs, other names."""
    out = {k: [] for k in ("id", "modes", "ra", "dec", "ogle4", "ogle3", "ogle2", "other")}
    for line in (C.OGLE_DIR / gal / "ident.dat").read_text().splitlines():
        if not line.strip():
            continue
        out["id"].append(line[0:17].strip())
        out["modes"].append(line[18:26].strip())
        ra_h = int(line[28:30]) + int(line[31:33]) / 60 + float(line[34:39]) / 3600
        out["ra"].append(15.0 * ra_h)
        out["dec"].append(_dms(line[40], line[41:43], line[44:46], line[47:51]))
        out["ogle4"].append(line[53:69].strip())
        out["ogle3"].append(line[70:85].strip())
        out["ogle2"].append(line[86:101].strip())
        out["other"].append(line[102:].strip())
    return {k: np.array(v) for k, v in out.items()}


def harvard_numbers(other: str) -> list[int]:
    """Harvard Variable ('HV') numbers in an OGLE 'other designation' field."""
    return [int(tok[2:]) for tok in other.split() if tok.startswith("HV") and tok[2:].isdigit()]


@lru_cache(maxsize=None)
def read_remarks(gal: str) -> dict:
    path = C.OGLE_DIR / gal / "remarks.txt"
    out = {}
    if path.exists():
        for line in path.read_text(errors="replace").splitlines():
            if line.strip():
                out[line[:17].strip()] = line[17:].strip()
    return out


@lru_cache(maxsize=None)
def catalog(gal: str) -> dict:
    """One row per Cepheid: modes, position, mean V and I, periods per mode, Fourier terms.

    ``mode`` is OGLE's mode string with slashes removed ('F', '1O', 'F1O', ...).
    ``P_F``, ``P_1O``, ``P_2O``, ``P_3O`` hold the period of each excited mode (NaN
    if absent). ``P`` is the dominant (longest) period, and the amplitude and
    Fourier parameters (R21, phi21, R31, phi31) belong to it.
    """
    ident = read_ident(gal)
    index = {sid: k for k, sid in enumerate(ident["id"])}
    n = len(index)
    cols = {name: np.full(n, np.nan) for name in
            ("I", "V", "P", "eP", "T0", "ampI", "R21", "phi21", "R31", "phi31",
             "P_F", "P_1O", "P_2O", "P_3O", "eP_F", "eP_1O", "T0_F", "T0_1O")}
    mode = np.array([""] * n, dtype=object)
    for fname, modes in MODE_FILES.items():
        path = C.OGLE_DIR / gal / fname
        if not path.exists():
            continue
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            k = index[line[0:17].strip()]
            mode[k] = "".join(modes)
            cols["I"][k] = _num(line[19:25])
            cols["V"][k] = _num(line[26:32])
            for j, (m, block) in enumerate(zip(modes, (_BLOCK, _BLOCK2, _BLOCK3))):
                vals = {key: _num(line[a:b]) for key, (a, b) in block.items()}
                cols[f"P_{m}"][k] = vals["P"]
                if m in ("F", "1O"):
                    cols[f"eP_{m}"][k] = vals["eP"]
                    cols[f"T0_{m}"][k] = vals["T0"]
                if j == 0:
                    for key in ("P", "eP", "T0", "ampI", "R21", "phi21", "R31", "phi31"):
                        cols[key][k] = vals[key]
    if np.any(mode == ""):
        missing = ident["id"][mode == ""]
        raise ValueError(f"{gal}: {len(missing)} stars in ident.dat without parameters, e.g. {missing[:3]}")
    remarks = read_remarks(gal)
    out = {"id": ident["id"], "gal": np.array([gal] * n), "mode": mode.astype(str),
           "ra": ident["ra"], "dec": ident["dec"], "other": ident["other"],
           "remark": np.array([remarks.get(s, "") for s in ident["id"]])}
    out.update(cols)
    out["logP"] = np.log10(out["P"])
    out["W_I"] = out["I"] - C.R_WI * (out["V"] - out["I"])
    return out


def read_lightcurve(gal: str, star_id: str, band: str = "I"):
    """(HJD - 2450000, mag, err) arrays, or None if OGLE has no light curve."""
    path = C.OGLE_DIR / gal / "phot" / band / f"{star_id}.dat"
    if not path.exists():
        return None
    a = np.loadtxt(path, ndmin=2)
    if a.size == 0:
        return None
    return a[:, 0], a[:, 1], a[:, 2]


@lru_cache(maxsize=None)
def galactic_list() -> dict:
    """OGLE's list of known Milky Way classical Cepheids (Pietrukowicz et al. 2021, 2026 update)."""
    out = {k: [] for k in ("source", "name", "ra", "dec", "l", "b", "mode", "P_F", "P_1O", "P_2O",
                           "P_3O", "bands", "ogle", "source_id", "G")}
    for line in (C.OGLE_DIR / C.OGLE_GALACTIC_LIST).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        out["source"].append(t[0])
        out["name"].append(re.sub("_+", " ", t[1]).strip())
        out["ra"].append(float(t[4]))
        out["dec"].append(float(t[5]))
        out["l"].append(float(t[6]))
        out["b"].append(float(t[7]))
        out["mode"].append(t[8])
        for key, tok in zip(("P_F", "P_1O", "P_2O", "P_3O"), t[9:13]):
            v = float(tok)
            out[key].append(v if v > 0 else np.nan)
        out["bands"].append(t[13])
        out["ogle"].append(t[14])
        out["source_id"].append(int(t[-2]) if t[-2].isdigit() else 0)
        out["G"].append(float(t[-1]) if t[-1] != "-" else np.nan)
    return {k: np.array(v) for k, v in out.items()}
