"""Period-eccentricity diagrams and tidal circularization.

Two stars on a close orbit raise tides on each other. Tidal friction drains
orbital energy while conserving angular momentum, and the orbit that keeps the
most angular momentum for its energy is a circle: close binaries end up on
circular orbits. The longer the stars have had, the wider the orbits that have
been circularized, so the *circularization period* is a clock (Meibom &
Mathieu 2005).

We measure it two ways on each sample:

* the Meibom & Mathieu (2005) recipe: fit e(P) = alpha (1 - exp(beta (P' - P)))
  for P > P' (and e = 0 below), with their alpha and beta, and quote the period
  where the fitted curve reaches e = 0.01;
* a model-light check: the period at which half the orbits are consistent
  with circular (Lucy & Sweeney 1971: an eccentricity is significant at 5% only
  if e > 2.45 sigma_e).

Gaia's SB1 eccentricities, like all fitted eccentricities, are biased upward
for near-circular orbits measured at low signal-to-noise: noise can only make
a circle look elongated. The circular-fraction check absorbs that bias by
using each orbit's own uncertainty.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy import optimize

import binaries as B
import samples
import stellar

#: Solar motion relative to the local standard of rest (Schoenrich, Binney &
#: Dehnen 2010, MNRAS 403, 1829), km/s.
SOLAR_UVW = (11.1, 12.24, 7.25)

#: Kinematic populations by total velocity relative to the LSR (the Toomre-
#: diagram convention of Nissen 2004 and Venn et al. 2004): thin disc stars are
#: mostly younger than ~8 Gyr, thick-disc and halo stars older than ~10 Gyr.
POPULATIONS = {"thin disc": (0, 50), "thick disc": (70, 180), "halo": (200, 1000)}

SB_TYPES = ("SB1", "SB1C", "SB2", "SB2C", "AstroSpectroSB1")


def _f(col):
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


def gaia_spectroscopic_orbits() -> Table:
    """DR3 spectroscopic orbits with absolute magnitudes and kinematics."""
    o = samples.gaia_table("nss_orbits")
    t = o[np.isin(np.char.strip(o["nss_solution_type"].astype(str)), SB_TYPES)].copy()
    t["nss_solution_type"] = np.char.strip(t["nss_solution_type"].astype(str))
    plx = _f(t["parallax"])
    nss_plx = _f(t["nss_parallax"])
    t["plx_best"] = np.where(np.isfinite(nss_plx), nss_plx, plx)
    t["parallax"] = t["plx_best"]
    samples.add_photometry(t)
    t["teff"] = stellar.teff_from_bprp(t["bp_rp"])
    t["circular_model"] = np.isin(t["nss_solution_type"], ("SB1C", "SB2C"))
    ecc = _f(t["eccentricity"])
    t["ecc"] = np.where(t["circular_model"], 0.0, ecc)
    t["ecc_err"] = np.where(t["circular_model"], 0.0, _f(t["eccentricity_error"]))
    add_kinematics(t)
    return t


def add_kinematics(t: Table) -> None:
    """Galactic space velocities (U, V, W) relative to the LSR, km/s."""
    rv = _f(t["center_of_mass_velocity"]) if "center_of_mass_velocity" in t.colnames else _f(t["radial_velocity"])
    plx = _f(t["parallax"])
    ok = np.isfinite(rv) & (plx > 0) & np.isfinite(_f(t["pmra"])) & np.isfinite(_f(t["pmdec"]))
    n = len(t)
    uvw = np.full((3, n), np.nan)
    if ok.any():
        c = SkyCoord(ra=_f(t["ra"])[ok] * u.deg, dec=_f(t["dec"])[ok] * u.deg,
                     distance=(1000.0 / plx[ok]) * u.pc, pm_ra_cosdec=_f(t["pmra"])[ok] * u.mas / u.yr,
                     pm_dec=_f(t["pmdec"])[ok] * u.mas / u.yr, radial_velocity=rv[ok] * u.km / u.s,
                     frame="icrs")
        g = c.galactic.cartesian.differentials["s"]
        uvw[:, ok] = np.vstack([g.d_x.to_value(u.km / u.s), g.d_y.to_value(u.km / u.s),
                                g.d_z.to_value(u.km / u.s)])
    t["U"] = uvw[0] + SOLAR_UVW[0]
    t["V"] = uvw[1] + SOLAR_UVW[1]
    t["W"] = uvw[2] + SOLAR_UVW[2]
    t["v_tot"] = np.sqrt(t["U"] ** 2 + t["V"] ** 2 + t["W"] ** 2)


def population(t: Table) -> np.ndarray:
    out = np.full(len(t), "", dtype="U10")
    for name, (lo, hi) in POPULATIONS.items():
        out[(t["v_tot"] >= lo) & (t["v_tot"] < hi)] = name
    return out


def consistent_with_circular(ecc, ecc_err, *, k: float = 2.45, floor: float = 0.0):
    """Lucy & Sweeney (1971): e is significant at 5% only if e > 2.45 sigma_e."""
    ecc = np.asarray(ecc, dtype=float)
    err = np.asarray(ecc_err, dtype=float)
    return (ecc <= np.maximum(k * np.nan_to_num(err, nan=0.0), floor))


# ---------------------------------------------------------------------------
# Circularization period estimators
# ---------------------------------------------------------------------------

def fit_meibom_mathieu(period, ecc, *, alpha=0.35, beta=0.14, p_max=None, p_grid=None) -> dict:
    """Best P' by least absolute deviation of e from the circularization
    function (alpha, beta fixed), and the circularization period e(P) = 0.01."""
    period = np.asarray(period, dtype=float)
    ecc = np.asarray(ecc, dtype=float)
    ok = np.isfinite(period) & np.isfinite(ecc)
    if p_max is not None:
        ok &= period < p_max
    p, e = period[ok], ecc[ok]
    if p_grid is None:
        p_grid = np.arange(0.5, 40.0, 0.01)
    cost = np.array([np.sum(np.abs(e - B.meibom_mathieu_e(p, pp, alpha, beta))) for pp in p_grid])
    k = int(np.argmin(cost))
    p_prime = float(p_grid[k])
    return {"p_prime": p_prime, "p_circ": float(B.meibom_mathieu_pcirc(p_prime, alpha, beta)),
            "n": int(ok.sum()), "cost": cost, "p_grid": p_grid}


def circular_fraction_curve(period, is_circular, *, edges=None) -> dict:
    """Fraction of orbits consistent with circular in bins of log P, and the
    period where a logistic curve through them crosses one half."""
    logp = np.log10(np.asarray(period, dtype=float))
    circ = np.asarray(is_circular, dtype=bool)
    if edges is None:
        edges = np.arange(-0.4, 2.01, 0.1)
    idx = np.digitize(logp, edges) - 1
    nb = len(edges) - 1
    n = np.bincount(idx[(idx >= 0) & (idx < nb)], minlength=nb)
    k = np.bincount(idx[(idx >= 0) & (idx < nb)], weights=circ[(idx >= 0) & (idx < nb)], minlength=nb)

    # logistic in log P with a floor (spurious/near-circular long-period orbits)
    def nll(p):
        x50, width, floor, top = p
        if not (0 < width < 2 and 0 <= floor < 0.5 and 0.5 < top <= 1):
            return 1e12
        f = floor + (top - floor) / (1 + np.exp((logp - x50) / width))
        f = np.clip(f, 1e-9, 1 - 1e-9)
        return -np.sum(np.where(circ, np.log(f), np.log(1 - f)))

    sel = np.isfinite(logp) & (logp > edges[0]) & (logp < edges[-1])
    logp_all, circ_all = logp, circ
    logp, circ = logp_all[sel], circ_all[sel]
    best = None
    for x0 in (0.6, 0.8, 1.0, 1.2):
        r = optimize.minimize(nll, [x0, 0.1, 0.1, 0.95], method="Nelder-Mead",
                              options={"maxiter": 5000, "xatol": 1e-5, "fatol": 1e-6})
        if best is None or r.fun < best.fun:
            best = r
    x50, width, floor, top = best.x
    # P where the circular fraction is halfway between its plateau and floor
    return {"edges": edges, "n": n, "k": k, "frac": np.where(n > 0, k / np.maximum(n, 1), np.nan),
            "logp50": float(x50), "p50": float(10 ** x50), "width": float(width), "floor": float(floor),
            "top": float(top), "n_used": int(sel.sum())}


def bootstrap(fn, arrays, n_boot=300, seed=1) -> np.ndarray:
    """Resample objects with replacement and re-run an estimator."""
    rng = np.random.default_rng(seed)
    n = len(arrays[0])
    out = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        out.append(fn(*[a[i] for a in arrays]))
    return np.array(out)


def fit_envelope(period, ecc, *, p_grid=None, tau_grid=None, delta=0.09):
    """Bashi, Mazeh & Faigler (2023)'s estimator: the upper envelope
    e_max(P) = 1 - (P_cut / P)^(1/tau) (zero below P_cut) of significantly
    eccentric orbits, with a soft (Fermi-function) edge of width delta; the
    density below the envelope is uniform in e. P_cut is the base of the
    envelope: the shortest period at which eccentric orbits are allowed.
    Maximum likelihood on a grid (their delta = 0.09 is fixed)."""
    period = np.asarray(period, dtype=float)
    ecc = np.asarray(ecc, dtype=float)
    ok = np.isfinite(period) & np.isfinite(ecc)
    period, ecc = period[ok], ecc[ok]
    if len(period) < 10:
        return {"p_cut": np.nan, "tau": np.nan, "n": int(len(period))}
    if p_grid is None:
        p_grid = np.arange(1.0, 12.0, 0.05)
    if tau_grid is None:
        tau_grid = np.arange(1.0, 3.01, 0.1)
    best = (-np.inf, np.nan, np.nan)
    for tau in tau_grid:
        for pc in p_grid:
            env = np.where(period > pc, 1 - (pc / period) ** (1 / tau), 0.0)
            # log of the Fermi edge, and its integral over e in [0, 1] in closed form:
            # int_0^1 de / (1 + exp((e - a)/d)) = 1 - d ln(1 + e^((1-a)/d)) + d ln(1 + e^(-a/d))
            log_f = -np.logaddexp(0.0, (ecc - env) / delta)
            norm = 1.0 - delta * np.logaddexp(0.0, (1 - env) / delta) + delta * np.logaddexp(0.0, -env / delta)
            ll = float(np.sum(log_f - np.log(norm)))
            if ll > best[0]:
                best = (ll, pc, tau)
    return {"p_cut": float(best[1]), "tau": float(best[2]), "n": int(len(period)), "lnl": best[0]}
