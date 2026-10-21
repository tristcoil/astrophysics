"""Weighing the planets by the star's wobble: the Hipparcos-Gaia accelerations.

The star circles the system's barycentre in the opposite sense to its planets:

    r_star(t) = - sum_k (m_k / M_total) X_k(t)            (exact, X_k relative to the star)

The HGCA (Brandt 2021, EDR3 edition) gives three proper motions of HR 8799: the
Hipparcos one (1989.9-1993.2), the Gaia EDR3 one (2014.6-2017.4) and the
"long-term" one from the change of position between the two missions (~1991.3
to ~2015.9). The barycentre's own motion cancels in the differences

    Delta_G = mu_Gaia - mu_HG,    Delta_H = mu_Hip - mu_HG,

which depend only on the reflex orbit: linear in the planet masses, with
coefficients set by the orbits. Brandt et al. (2021) used this to weigh
HR 8799 e (9.6 +1.9/-1.8 MJup, given mass ratios between the planets).

Two tests here: (1) the mass scale implied by the measured Delta_G and
Delta_H for given orbits and mass ratios; (2) whether the direction of the
measured acceleration agrees with the four planets at all - Lagrange et al.
(2026) argue that a fifth, inner companion is needed.
"""

from __future__ import annotations

import json

import numpy as np
from astropy.time import Time

import config as C
import kepler as K

#: mission windows (decimal year): Hipparcos 1989.85-1993.21, Gaia EDR3 2014.57-2017.41
HIP_WINDOW = (1989.85, 1993.21)
GAIA_WINDOW = (2014.57, 2017.41)


def load_hgca() -> dict:
    return json.loads((C.DERIVED / "hgca_hr8799.json").read_text())


def measured() -> tuple[np.ndarray, np.ndarray]:
    """(Delta_G_ra, Delta_G_dec, Delta_H_ra, Delta_H_dec) in mas/yr and covariance.

    The long-term (HG) proper motion is common to both differences, so their
    errors are correlated; within each proper motion we use the HGCA's
    ra-dec correlation coefficients.
    """
    h = load_hgca()

    def cov(pref):
        sa, sd, rho = h[f"pmra_{pref}_error"], h[f"pmdec_{pref}_error"], h[f"pmra_pmdec_{pref}"]
        return np.array([[sa * sa, rho * sa * sd], [rho * sa * sd, sd * sd]])

    g = np.array([h["pmra_gaia"], h["pmdec_gaia"]])
    hg = np.array([h["pmra_hg"], h["pmdec_hg"]])
    hip = np.array([h["pmra_hip"], h["pmdec_hip"]])
    d = np.concatenate([g - hg, hip - hg])
    Cg, Chg, Ch = cov("gaia"), cov("hg"), cov("hip")
    big = np.zeros((4, 4))
    big[:2, :2] = Cg + Chg
    big[2:, 2:] = Ch + Chg
    big[:2, 2:] = Chg
    big[2:, :2] = Chg
    return d, big


def _mjd(year):
    return Time(year, format="decimalyear").mjd


def reflex_design(elements: dict, plx, mstar, masses=None, n_window: int = 9) -> np.ndarray:
    """Per planet k, the change in (Delta_G, Delta_H) per Msun of planet mass.

    elements: planet -> (a, e, inc, aop, pan, tau) arrays of shape (n,) (Jacobi);
    returns (n, 4 planets, 4) in mas/yr per Msun, before dividing by M_total.
    Proper motions are modelled as the mean velocity over each mission (a
    straight-line fit to positions sampled through the window), the HG one as
    the position change between the two missions' central epochs.
    """
    h = load_hgca()
    tg = np.linspace(*GAIA_WINDOW, n_window)
    th = np.linspace(*HIP_WINDOW, n_window)
    ep = {"g": tg, "h": th}
    # central epochs of the HG difference (per coordinate in the HGCA)
    t_hg = {"ra": (h["epoch_ra_hip"], h["epoch_ra_gaia"]), "dec": (h["epoch_dec_hip"], h["epoch_dec_gaia"])}
    times = np.concatenate([tg, th, [t_hg["ra"][0], t_hg["ra"][1], t_hg["dec"][0], t_hg["dec"][1]]])
    mjd = _mjd(times)[None, :]
    el = {p: tuple(np.asarray(v, float)[:, None] for v in elements[p]) for p in C.INNER_TO_OUTER}
    sky = K.system_offsets(mjd, el, np.asarray(plx, float)[:, None], np.asarray(mstar, float)[:, None],
                           masses)
    n = len(np.atleast_1d(plx))
    out = np.zeros((n, 4, 4))
    for k, p in enumerate(C.INNER_TO_OUTER):
        ra, dec = sky[p]
        # star reflex per unit (m_k / M_total): minus the planet's offset
        sra, sdec = -ra, -dec
        ng = len(tg)
        fit = lambda t, y: np.polyfit(t - t.mean(), y.T, 1)[0]   # slope per sample
        mu_g = np.stack([fit(tg, sra[:, :ng]), fit(tg, sdec[:, :ng])], axis=1)
        mu_h = np.stack([fit(th, sra[:, ng:2 * ng]), fit(th, sdec[:, ng:2 * ng])], axis=1)
        i0 = 2 * ng
        mu_hg = np.stack([(sra[:, i0 + 1] - sra[:, i0]) / (t_hg["ra"][1] - t_hg["ra"][0]),
                          (sdec[:, i0 + 3] - sdec[:, i0 + 2]) / (t_hg["dec"][1] - t_hg["dec"][0])], axis=1)
        out[:, k, :2] = mu_g - mu_hg
        out[:, k, 2:] = mu_h - mu_hg
    return out


def mass_scale(design: np.ndarray, mtot: np.ndarray, ratios: np.ndarray, data=None) -> dict:
    """Best mass scale s (MJup) with masses = s * ratios, per orbit sample.

    design: (n, 4, 4) from reflex_design; ratios: (4,) or (n, 4) in the order
    e, d, c, b. Returns per-sample best s, its error, and the chi^2 of the fit
    (2 measured vectors, 1 free scale: 3 degrees of freedom).
    """
    d, cov = measured() if data is None else data
    W = np.linalg.inv(cov)
    r = np.broadcast_to(ratios, (design.shape[0], 4))
    G = np.einsum("nkj,nk->nj", design, r) * C.MJUP_MSUN / mtot[:, None]   # per MJup of scale
    a = np.einsum("ni,ij,nj->n", G, W, G)
    b = np.einsum("ni,ij,j->n", G, W, d)
    s = b / a
    res = d[None, :] - s[:, None] * G
    chi2 = np.einsum("ni,ij,nj->n", res, W, res)
    return {"scale": s, "scale_err": 1 / np.sqrt(a), "chi2": chi2, "G": G}


def predicted(design: np.ndarray, mtot: np.ndarray, masses_mjup: np.ndarray) -> np.ndarray:
    """Predicted (Delta_G, Delta_H) in mas/yr for given masses (n, 4) in MJup."""
    return np.einsum("nkj,nk->nj", design, masses_mjup * C.MJUP_MSUN) / mtot[:, None]
