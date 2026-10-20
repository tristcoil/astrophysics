"""Rebuild Raghavan et al. (2010)'s period distribution from their own table.

Their Table 13 lists, for each of the 454 Sun-like stars within 25 pc, every
companion with its period, or its separation, and how it was found. Their
Fig. 13 is a histogram of log P for the 259 confirmed stellar and brown-dwarf
companions (planets, also in the table, excluded), with a Gaussian fit:
mean log P = 5.03, sigma = 2.28 (P in days). Periods of pairs without orbits
were estimated as they describe (Sec. 5.3.3):

* resolved pairs: log a" = log rho" + 0.13 (Duquennoy & Mayor 1991), then
  Kepler's third law with the stars' masses;
* unresolved companions with no orbit: 30-200 yr for radial-velocity
  variables, 10-25 yr for proper-motion accelerations (we take the
  logarithmic middle of each range).

Rebuilding it lets us check that we read their result correctly, and gives us
their companions per primary per decade of period -- the quantity our Gaia
census measures.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

import binaries as B
import config as C
import literature as L

YEAR = 365.25


def _s(col):
    return np.char.strip(np.ma.filled(np.ma.asarray(col).astype(str), ""))


def _f(col):
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


def companion_periods() -> Table:
    t = Table.read(C.VIZIER_DIR / "raghavan10_companions.fits")
    stars = Table.read(C.VIZIER_DIR / "raghavan10_stars.fits")
    st, comp = _s(t["St"]), _s(t["Comp"])
    is_comp = (st == "Y") & (comp != "")
    per, unit = _f(t["Per"]), _s(t["x_Per"])
    asep, lsep = _f(t["Asep"]), _f(t["Lsep"])
    sb, ot = _s(t["SB"]), _s(t["OT"])
    hip = np.ma.filled(np.ma.asarray(t["HIP"]), -1)
    mass_by_hip = {int(h): float(m) for h, m in zip(np.ma.filled(stars["HIP"], -1), _f(stars["Mass"]))
                   if np.isfinite(m)}

    p_days = np.where(unit == "y", per * YEAR, np.where(unit == "h", per / 24.0, per))
    how = np.where(np.isfinite(p_days), "orbit", "")
    m1 = np.array([mass_by_hip.get(int(h), 1.0) for h in hip])
    m_tot = m1 * 1.6  # companion mass unknown for most wide pairs; P ~ M^-1/2 makes this a 0.1 dex effect
    # resolved pairs without an orbit: a = 10^0.13 x projected separation
    sep_ok = ~np.isfinite(p_days) & np.isfinite(lsep) & np.isfinite(asep)
    p_days = np.where(sep_ok, B.period_days_from_a(L.RAGHAVAN_A_OVER_S * lsep, m_tot), p_days)
    how = np.where(sep_ok, "separation", how)
    # Lsep without Asep is a semi-major axis from an orbit (ReadMe note 3)
    a_ok = ~np.isfinite(p_days) & np.isfinite(lsep) & ~np.isfinite(asep)
    p_days = np.where(a_ok, B.period_days_from_a(lsep, m_tot), p_days)
    how = np.where(a_ok, "semi-major axis", how)
    rv = ~np.isfinite(p_days) & (sb == "V")
    p_days = np.where(rv, np.sqrt(30 * 200) * YEAR, p_days)
    how = np.where(rv, "RV variable (30-200 yr)", how)
    acc = ~np.isfinite(p_days) & (ot == "M")
    p_days = np.where(acc, np.sqrt(10 * 25) * YEAR, p_days)
    how = np.where(acc, "acceleration (10-25 yr)", how)

    out = Table({"hip": hip, "comp": comp, "period_days": p_days, "how": how, "sb": sb, "ot": ot,
                 "lsep_au": lsep})
    return out[is_comp]


def refit() -> dict:
    t = companion_periods()
    lp = np.log10(np.asarray(t["period_days"]))
    ok = np.isfinite(lp)
    lp = lp[ok]
    # unbinned maximum likelihood = sample mean and standard deviation
    mean, sd = float(np.mean(lp)), float(np.std(lp))
    # least-squares Gaussian to a 1-dex histogram (the binning of their Fig. 13)
    edges = np.arange(-1.0, 11.01, 1.0)
    h, _ = np.histogram(lp, edges)
    c = 0.5 * (edges[1:] + edges[:-1])
    from scipy.optimize import curve_fit
    g = lambda x, a, mu, s: a * np.exp(-0.5 * ((x - mu) / s) ** 2)
    (a, mu, s), _ = curve_fit(g, c, h, p0=[h.max(), 5, 2.3])
    return {"n_companions": int(len(t)), "n_with_period": int(ok.sum()),
            "by_method": {k: int(np.sum(t["how"] == k)) for k in np.unique(t["how"])},
            "ml_mean": mean, "ml_sigma": sd, "hist_mean": float(mu), "hist_sigma": float(abs(s)),
            "hist_edges": edges, "hist": h, "logp": lp,
            # companions per primary per decade, the quantity Gaia measures
            "per_primary_per_dex": h / L.RAGHAVAN_N_PRIMARIES / np.diff(edges)}
