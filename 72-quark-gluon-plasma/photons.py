"""The plasma starts hotter: thermal photons.

Photons interact only electromagnetically, so once made they leave the
fireball unscattered and carry the temperature of the moment they were
emitted. ALICE measured the direct-photon spectrum (photons not from hadron
decays) in central Pb-Pb at 2.76 TeV; between 0.9 and 2.1 GeV/c, where
prompt photons from the initial hard scatterings are small, it falls as
exp(-pT / T_eff). T_eff is a space-time average over the whole evolution and
is blue-shifted by the fireball's expansion, so it is not the initial
temperature; hydrodynamic models that reproduce it start at 300-600 MeV.

We fit the published 0-20% spectrum (HEPData ins1394677 Table 1) with the
statistical errors, and get the systematic error by tilting all points by
their (correlated) systematic errors across the fit range.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

import config as cfg
import datasets


def spectrum(table: str = "Table1.yaml") -> dict:
    t = datasets.hepdata_table("alice_photons", table)
    lo = np.array([x["low"] for x in t["independent_variables"][0]["values"]], float)
    hi = np.array([x["high"] for x in t["independent_variables"][0]["values"]], float)
    vals, stat, sys = [], [], []
    for v in t["dependent_variables"][0]["values"]:
        val = v["value"]
        if isinstance(val, str):                    # upper limits ("<0.645")
            vals.append(np.nan), stat.append(np.nan), sys.append(np.nan)
            continue
        s, y = datasets._errors(float(val), v.get("errors"))
        vals.append(float(val)), stat.append(s), sys.append(y)
    return {"lo": lo, "hi": hi, "pt": 0.5 * (lo + hi), "y": np.array(vals),
            "stat": np.array(stat), "sys": np.array(sys)}


def _fit(pt, y, s):
    f = lambda x, a, T: a * np.exp(-x / T)          # noqa: E731
    p, c = curve_fit(f, pt, y, p0=(y[0] * np.exp(pt[0] / 0.3), 0.3), sigma=s,
                     absolute_sigma=True, maxfev=20000)
    chi2 = float(np.sum(((y - f(pt, *p)) / s) ** 2))
    return p, np.sqrt(np.diag(c)), chi2


def inverse_slope(fit_range: tuple = cfg.ALICE_PHOTONS_2016["fit_range"]) -> dict:
    sp = spectrum()
    m = (sp["lo"] >= fit_range[0] - 1e-9) & (sp["hi"] <= fit_range[1] + 1e-9) & np.isfinite(sp["y"])
    pt, y, st, sy = sp["pt"][m], sp["y"][m], sp["stat"][m], sp["sys"][m]
    p, e, chi2 = _fit(pt, y, st)
    # Systematic error: refit with statistical and systematic errors in
    # quadrature (treated as uncorrelated); the extra width is the systematic.
    # A full +-1 sigma linear tilt of the spectrum would give +-90 MeV, an
    # extreme ALICE's partially correlated type-B errors do not imply.
    pt_all, tot = _fit(pt, y, np.hypot(st, sy))[0:2]
    syst = float(np.sqrt(max(tot[1] ** 2 - e[1] ** 2, 0.0)))
    return {"T_eff_MeV": p[1] * 1e3, "stat_MeV": e[1] * 1e3, "syst_MeV": syst * 1e3,
            "chi2": chi2, "ndf": int(m.sum()) - 2, "A": p[0], "n_points": int(m.sum()),
            "T_eff_total_errors_MeV": pt_all[1] * 1e3, "err_total_errors_MeV": tot[1] * 1e3,
            "pt": pt.tolist(), "y": y.tolist(), "stat": st.tolist(), "sys": sy.tolist(),
            "published": cfg.ALICE_PHOTONS_2016}
