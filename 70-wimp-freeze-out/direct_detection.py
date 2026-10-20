"""Direct detection: WIMPs bumping into xenon nuclei underground.

If the Milky Way's dark-matter halo is made of WIMPs, a few of them should occasionally scatter
off atomic nuclei in a detector. The leading experiments fill a time-projection chamber with tonnes
of liquid xenon and look for single nuclear recoils of a few keV. None has been seen. The limits
below are 90% CL upper limits on the spin-independent WIMP-nucleon cross section, all computed by
the collaborations for the standard halo (local density 0.3 GeV/cm^3, Maxwellian velocities):

* LZ 2024 (4.2 tonne-years; PRL 135, 011802), HEPData record ins2841863;
* XENONnT 2025 (SR0+SR1, 3.1 tonne-years; arXiv:2502.18005), the collaboration's data release;
* earlier limits from C. O'Hare's compilation (LUX 2017, PandaX-II 2017, XENON1T 2018, ...);
* the xenon "neutrino floor" (O'Hare 2021): where solar and atmospheric neutrinos scattering
  coherently on nuclei start to look like WIMPs.

On this plane we put one thermal relic: the Higgs-portal scalar singlet (higgs_portal.py), whose
freeze-out fixes its coupling and therefore its scattering rate.
"""

from __future__ import annotations

import csv
import json
import math

import numpy as np

import config as C

DD_DIR = C.RAW_DIR / "dd"


def _blocks(path):
    """HEPData CSV with several tables separated by blank lines: {column name: (x, y)}."""
    out, cur, header = {}, [], None
    for line in path.read_text().splitlines() + [""]:
        if line.startswith("#"):
            continue
        if not line.strip():
            if header and cur:
                a = np.array(cur, dtype=float)
                out[header[1]] = (a[:, 0], a[:, 1])
            cur, header = [], None
            continue
        if header is None:
            header = line.split(",")
            continue
        cur.append([float(v) for v in line.split(",")])
    return out


def lz2024() -> dict:
    b = _blocks(DD_DIR / "lz2024_si.csv")
    return {"limit": b["limit"], "median": b["median"], "minus1": b["-1sigma"], "plus1": b["1sigma"],
            "minus2": b["-2sigma"], "plus2": b["2sigma"]}


def xenonnt2025() -> dict:
    with open(DD_DIR / "xenonnt2025_si.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    m = np.array([float(r["mass"]) for r in rows])
    get = (lambda k: np.array([float(r[k]) for r in rows]))
    return {"limit": (m, get("upper_limit")), "median": (m, get("sensitivity_0")),
            "minus1": (m, get("sensitivity_-1")), "plus1": (m, get("sensitivity_1"))}


OHARE = {"LUX 2017": "ohare_LUX.txt", "PandaX-II 2017": "ohare_PandaX.txt",
         "XENON1T 2018": "ohare_XENON1T-2018.txt", "DarkSide-50": "ohare_DarkSide.txt",
         "CRESST": "ohare_CRESST.txt", "CDMSlite": "ohare_CDMSLite.txt"}


def ohare(name: str):
    d = np.loadtxt(DD_DIR / OHARE[name], comments="#")
    o = np.argsort(d[:, 0])
    return d[o, 0], d[o, 1]


def dama_regions():
    """DAMA/LIBRA's annual modulation interpreted as WIMPs (Savage et al. 2009): two islands."""
    return [np.loadtxt(DD_DIR / f, comments="#") for f in ("ohare_DAMA1.txt", "ohare_DAMA2.txt")]


def neutrino_floor():
    d = np.loadtxt(DD_DIR / "ohare_Xe_SI_floor.txt", comments="#")
    return d[:, 0], d[:, 1]


def limit_at(m, curve) -> np.ndarray:
    """Log-log interpolation; above the table, extrapolate as sigma ~ m (rate ~ 1/m at high mass)."""
    x, y = curve
    m = np.atleast_1d(np.asarray(m, dtype=float))
    v = np.exp(np.interp(np.log(m), np.log(x), np.log(y), left=np.nan, right=np.nan))
    hi = m > x[-1]
    v[hi] = y[-1] * m[hi] / x[-1]
    return v


# --- the Higgs-portal line, cached ---------------------------------------------------------------

HP_CACHE = C.DERIVED_DIR / "higgs_portal_relic.json"
HP_MASSES = np.unique(np.concatenate([
    np.geomspace(20.0, 50.0, 8), np.linspace(52.0, 66.0, 29), np.geomspace(70.0, 300.0, 14),
    np.geomspace(320.0, 3.0e4, 26)]))


def higgs_portal_line(force: bool = False) -> dict:
    import higgs_portal as hp

    if HP_CACHE.exists() and not force:
        d = json.loads(HP_CACHE.read_text())
        if np.allclose(d["mass"], HP_MASSES):
            return d
    rows = [hp.relic_coupling(float(m)) for m in HP_MASSES]
    d = {"mass": [r["m"] for r in rows], "lambda": [r["lambda"] for r in rows],
         "omega_h2": [r["omega_h2"] for r in rows], "sigma_si": [r["sigma_si_cm2"] for r in rows],
         "lambda_inv_max": [hp.lambda_invisible_max(float(m)) for m in HP_MASSES],
         "f_N": hp.F_N, "f_N_err": hp.F_N_ERR, "gamma_h_GeV": hp.gamma_h()}
    HP_CACHE.parent.mkdir(parents=True, exist_ok=True)
    HP_CACHE.write_text(json.dumps(d, indent=1))
    return d


def allowed_windows(mass, sigma_pred, limit_curve, extra_ok=None) -> list:
    """Mass intervals where the prediction lies below the limit (and extra_ok holds)."""
    m = np.asarray(mass)
    ok = np.asarray(sigma_pred) < limit_at(m, limit_curve)
    if extra_ok is not None:
        ok &= np.asarray(extra_ok)
    wins, start = [], None
    for i, flag in enumerate(ok):
        if flag and start is None:
            start = i
        if (not flag or i == len(ok) - 1) and start is not None:
            end = i if flag else i - 1
            wins.append([float(m[start]), float(m[end])])
            start = None
    return wins
