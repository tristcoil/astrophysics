"""Henrietta Leavitt's 25 Cepheids (Harvard College Observatory Circular 173, 1912).

Table I of the circular, transcribed by us from the scan the NASA ADS serves
(https://articles.adsabs.harvard.edu/pdf/1912HarCi.173....1L, cached in
data/raw/literature/) and checked digit by digit against a 400-dpi rendering.
Columns: Harvard number (H.), photographic magnitude at maximum and minimum
"as read from the light curve" on the Harvard provisional scale, epoch in days
after JD 2,410,000 (a reference time within the first cycle; the text does not
say it is a maximum), period in days, and Leavitt's residuals of maximum (M)
and minimum (m) from the straight lines of her Figure 2. ``new`` marks the
Harvard numbers printed in italics, "when the period has not been published
hitherto": the nine stars missing from Table VI of her 1908 catalogue (H.A. 60,
No. 4), which held 16 periods although the circular says seventeen and "8
additional". An independent reading of the same scan (and Shapley 1918, ApJ 48,
89, Table VII) agrees with ours in every entry. HV 1374's printed residual of
the minimum (-0.3) is probably a sign misprint (+0.3 fits her line); we keep it
as printed.

The modern counterparts are found through OGLE's own cross-identifications:
OGLE's ident.dat lists the Harvard Variable ('HV') number of 1,260 SMC Cepheids.
"""

from __future__ import annotations

import numpy as np

import config as C
import ogle

#: (H., max, min, epoch [JD - 2410000], period [d], res M, res m, italic/new)
TABLE_I = (
    (1505, 14.8, 16.1, 0.02, 1.25336, -0.6, -0.5, False),
    (1436, 14.8, 16.4, 0.02, 1.6637, -0.3, +0.1, False),
    (1446, 14.8, 16.4, 1.38, 1.7620, -0.3, +0.1, False),
    (1506, 15.1, 16.3, 1.08, 1.87502, +0.1, +0.1, False),
    (1413, 14.7, 15.6, 0.35, 2.17352, -0.2, -0.5, True),
    (1460, 14.4, 15.7, 0.00, 2.913, -0.3, -0.1, True),
    (1422, 14.7, 15.9, 0.6, 3.501, +0.2, +0.2, True),
    (842, 14.6, 16.1, 2.61, 4.2897, +0.3, +0.6, False),
    (1425, 14.3, 15.3, 2.8, 4.547, 0.0, -0.1, False),
    (1742, 14.3, 15.5, 0.95, 4.9866, +0.1, +0.2, False),
    (1646, 14.4, 15.4, 4.30, 5.311, +0.3, +0.1, False),
    (1649, 14.3, 15.2, 5.05, 5.323, +0.2, -0.1, False),
    (1492, 13.8, 14.8, 0.6, 6.2926, -0.2, -0.4, True),
    (1400, 14.1, 14.8, 4.0, 6.650, +0.2, -0.3, False),
    (1355, 14.0, 14.8, 4.8, 7.483, +0.2, -0.2, True),
    (1374, 13.9, 15.2, 6.0, 8.397, +0.2, -0.3, False),
    (818, 13.6, 14.7, 4.0, 10.336, 0.0, 0.0, False),
    (1610, 13.4, 14.6, 11.0, 11.645, 0.0, 0.0, True),
    (1365, 13.8, 14.8, 9.6, 12.417, +0.4, +0.2, True),
    (1351, 13.4, 14.4, 4.0, 13.08, +0.1, -0.1, True),
    (827, 13.4, 14.3, 11.6, 13.47, +0.1, -0.2, False),
    (822, 13.0, 14.6, 13.0, 16.75, -0.1, +0.3, True),
    (823, 12.2, 14.1, 2.9, 31.94, -0.3, +0.4, False),
    (824, 11.4, 12.8, 4.0, 65.8, -0.4, -0.2, False),
    (821, 11.2, 12.1, 97.0, 127.0, -0.1, -0.4, False),
)

#: Leavitt's own statement of the slope: "The logarithm of the period
#: increases by about 0.48 for each increase of one magnitude in brightness."
LEAVITT_DLOGP_PER_MAG = 0.48

#: Decimal places Leavitt quoted for each period (text strings above lose
#: trailing zeros), used as her implied precision: +-1 in the last digit.
PERIOD_DECIMALS = {1505: 5, 1436: 4, 1446: 4, 1506: 5, 1413: 5, 1460: 3, 1422: 3, 842: 4, 1425: 3, 1742: 4,
                   1646: 3, 1649: 3, 1492: 4, 1400: 3, 1355: 3, 1374: 3, 818: 3, 1610: 3, 1365: 3, 1351: 2,
                   827: 2, 822: 2, 823: 2, 824: 1, 821: 1}


def table() -> dict:
    a = np.array([r[:7] for r in TABLE_I], float)
    out = {"H": a[:, 0].astype(int), "max": a[:, 1], "min": a[:, 2], "epoch": a[:, 3], "P": a[:, 4],
           "res_max": a[:, 5], "res_min": a[:, 6], "new": np.array([r[7] for r in TABLE_I])}
    out["e_P"] = np.array([10.0 ** -PERIOD_DECIMALS[h] for h in out["H"]])
    out["logP"] = np.log10(out["P"])
    return out


def cross_identify(smc: dict) -> dict:
    """OGLE counterpart of each of Leavitt's stars, via the HV numbers OGLE lists."""
    t = table()
    hv_index = {}
    for k, other in enumerate(smc["other"]):
        for h in ogle.harvard_numbers(other):
            hv_index.setdefault(h, []).append(k)
    idx = np.full(len(t["H"]), -1)
    n_match = np.zeros(len(t["H"]), int)
    for j, h in enumerate(t["H"]):
        ks = hv_index.get(int(h), [])
        n_match[j] = len(ks)
        if ks:
            idx[j] = ks[0]
    t["ogle_index"] = idx
    t["n_ogle_matches"] = n_match
    return t


def compare(smc: dict) -> dict:
    """Leavitt's numbers against OGLE's for the same stars, and her law against today's."""
    import plfit
    t = cross_identify(smc)
    k = t["ogle_index"]
    ok = k >= 0
    get = lambda key: np.where(ok, smc[key][np.maximum(k, 0)], np.nan)
    t["ogle_id"] = np.where(ok, smc["id"][np.maximum(k, 0)], "")
    t["ogle_mode"] = np.where(ok, smc["mode"][np.maximum(k, 0)], "")
    t["P_ogle"] = get("P")
    t["eP_ogle"] = get("eP")
    for key in ("V", "I", "lc_V_max", "lc_V_min", "lc_I_max", "lc_I_min", "lc_V_amp", "lc_I_amp", "lc_I_T_max", "W_I"):
        t[key] = get(key)
    t["dP"] = t["P_ogle"] - t["P"]
    t["dP_over_P"] = t["dP"] / t["P"]
    # a period change over a century, in seconds per year, between Leavitt's photographs
    # (mean epoch ~1900; the plates span the 1890s-1911) and OGLE-IV (2010-2020)
    years = 2015.0 - 1900.0
    t["pdot_s_per_yr"] = t["dP"] * 86400.0 / years
    t["pdot_err_s_per_yr"] = np.hypot(t["e_P"], t["eP_ogle"]) * 86400.0 / years
    # her law, refitted: magnitude vs log P for maxima and minima (no clipping; 25 stars)
    fit_max = plfit.fit_line(t["logP"], t["max"], nsigma=1e9)
    fit_min = plfit.fit_line(t["logP"], t["min"], nsigma=1e9)
    mean_pg = 0.5 * (t["max"] + t["min"])
    fit_mean = plfit.fit_line(t["logP"], mean_pg, nsigma=1e9)
    # the same stars today, in V (OGLE light-curve extrema from our Fourier fits)
    okv = ok & np.isfinite(t["lc_V_max"])
    fit_vmax = plfit.fit_line(t["logP"], t["lc_V_max"], keep=okv, nsigma=1e9)
    fit_v = plfit.fit_line(t["logP"], t["V"], keep=ok & np.isfinite(t["V"]), nsigma=1e9)
    # Leavitt's mean photographic magnitude against OGLE's V: offset and scale
    m = ok & np.isfinite(t["V"])
    scale = np.polyfit(t["V"][m], mean_pg[m], 1)
    return {"table": t, "fit_max": fit_max, "fit_min": fit_min, "fit_mean": fit_mean, "fit_vmax": fit_vmax,
            "fit_v": fit_v, "pg_vs_V": {"slope": float(scale[0]), "intercept": float(scale[1]),
                                        "rms": float(np.std(mean_pg[m] - np.polyval(scale, t["V"][m])))},
            "n_matched": int(ok.sum()), "years": years}


#: Sampling frequencies that produce aliases in ground-based photographic
#: monitoring: the solar day, the year (observing seasons) and the synodic
#: month (plates are not taken near full moon), in cycles per day.
ALIASES = {"1 cycle/day": 1.0, "1 cycle/year": 1 / 365.2422, "1 cycle/lunar month": 1 / 29.5306,
           "2 cycles/year": 2 / 365.2422}


def aliases(t: dict, tol: float = 0.06) -> list:
    """Stars whose 1912 frequency differs from OGLE's by a sampling frequency (within ``tol``)."""
    out = []
    for j in range(len(t["H"])):
        if not np.isfinite(t["P_ogle"][j]):
            continue
        df = abs(1 / t["P"][j] - 1 / t["P_ogle"][j])
        # significant only if far larger than the quoted precision allows
        df_err = t["e_P"][j] / t["P"][j] ** 2
        if df < 5 * df_err:
            continue
        best = min(ALIASES.items(), key=lambda kv: abs(df / kv[1] - 1))
        frac = df / best[1] - 1
        out.append({"H": int(t["H"][j]), "P1912": float(t["P"][j]), "P_ogle": float(t["P_ogle"][j]),
                    "df": float(df), "alias": best[0] if abs(frac) < tol else "", "mismatch": float(frac)})
    return out
