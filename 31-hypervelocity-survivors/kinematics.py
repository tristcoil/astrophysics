"""From Gaia astrometry plus a radial velocity to a 3D velocity, a flight time
and a launch speed -- with honest distance uncertainties.

For most of these stars the proper motion (the motion across the sky) and the
radial velocity (along the line of sight, from the spectrum's Doppler shift)
are measured well. The distance is not: several parallaxes are consistent with
zero. Tangential speed is 4.74 x proper motion x distance, so the distance
uncertainty is the whole story, and the answer depends on what we assume about
distances before looking at the parallax (the *prior*). Following El-Badry et
al. (2023), section 5, we sample the posterior of

    distance d, proper motion (mu_a*, mu_d), radial velocity RV

given Gaia's parallax (corrected for its known zero-point offset, Lindegren et
al. 2021), the proper motions with their correlations, and the published RV,
under one of three distance priors (flat in absolute magnitude, the fiducial;
flat in distance; or an exponentially decreasing space density), with one
physical constraint: no launch faster than 3,000 km/s (two 1.3 Msun white
dwarfs in contact could not do better).

Because the parallax depends only on d, and the proper motions depend on d
only through their correlation with the parallax, the posterior factorises:
d is drawn from a 1D grid, then (mu, RV) from Gaussians conditioned on it, and
draws that break the 3,000 km/s limit are rejected. This gives exact,
independent samples (no Markov chain to converge).

Then, as in El-Badry et al., each sample's straight-line trajectory is traced
to the Galactic midplane (flight time t = z / v_z). The star was presumably
orbiting with the disc before the explosion, so the launch ("ejection") speed
is the present velocity minus the circular velocity at the crossing point.
"""

from __future__ import annotations

import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import CartesianDifferential, Galactocentric, ICRS, SkyCoord

from config import (A_G_PER_EBV, D_GRID_MAX, D_GRID_MIN, D_GRID_N, K_KMS, R0_KPC,
                    V_CIRC_FLAT, V_EJECTION_MAX, V_SUN_KMS, Z_SUN_KPC)

#: Casagrande & VandenBerg (2018): E(BP-RP) = (R_BP - R_RP) E(B-V).
E_BPRP_PER_EBV = 3.374 - 2.035
#: Length scale of the exponentially decreasing space density prior
#: (Astraatmadja & Bailer-Jones 2016), as used by El-Badry et al.
EDSD_L_KPC = 1.35


def galcen_frame() -> Galactocentric:
    return Galactocentric(galcen_distance=R0_KPC * u.kpc, z_sun=Z_SUN_KPC * u.kpc,
                          galcen_v_sun=CartesianDifferential(V_SUN_KMS * u.km / u.s))


def parallax_zero_point(row) -> float:
    """Lindegren et al. (2021) parallax zero point (mas) for one source; the
    *measured* parallax minus this is the corrected one. Falls back to the
    quasar median (-0.017 mas) outside the recipe's validity range."""
    try:
        from zero_point import zpt
        zpt.load_tables()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            z = zpt.get_zpt(np.array([float(row["phot_g_mean_mag"])]),
                            np.array([float(row["nu_eff_used_in_astrometry"])]),
                            np.array([float(row["pseudocolour"])]),
                            np.array([float(row["ecl_lat"])]),
                            np.array([int(row["astrometric_params_solved"])]), _warnings=False)
        z = float(z[0])
        return z if np.isfinite(z) else -0.017
    except Exception:  # noqa: BLE001
        return -0.017


# --------------------------------------------------------------------------
# Distance posterior
# --------------------------------------------------------------------------

def log_prior(d, prior: str, gauss=None):
    """ln p(d) up to a constant, d in kpc."""
    if prior == "flat_MG":          # flat in absolute magnitude == p(d) ~ 1/d
        return -np.log(d)
    if prior == "flat_d":
        return np.zeros_like(d)
    if prior == "edsd":
        return 2 * np.log(d) - d / EDSD_L_KPC
    if prior == "gauss":
        mu, sig = gauss
        return -0.5 * ((d - mu) / sig) ** 2
    raise ValueError(prior)


def distance_grid() -> np.ndarray:
    return np.geomspace(D_GRID_MIN, D_GRID_MAX, D_GRID_N)


def sample_distance(plx, plx_err, prior, n, rng, gauss=None) -> np.ndarray:
    """Draws of d (kpc) from p(d | parallax) on a log grid (with the Jacobian
    of the log spacing), by inverse-CDF sampling inside grid cells."""
    d = distance_grid()
    lnp = log_prior(d, prior, gauss) - 0.5 * ((plx - 1.0 / d) / plx_err) ** 2 + np.log(d)
    w = np.exp(lnp - lnp.max())
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]))])
    cdf /= cdf[-1]
    lnd = np.interp(rng.random(n), cdf, np.log(d))
    return np.exp(lnd)


# --------------------------------------------------------------------------
# Phase space and launch speed
# --------------------------------------------------------------------------

def _affine():
    """The ICRS -> Galactocentric transformation is a fixed rotation plus a
    shift (positions) and the Sun's velocity (velocities). Read it off
    astropy once, so the sampler can apply it as a matrix product."""
    if not hasattr(_affine, "cache"):
        pts = np.array([[1e-9, 1.0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1.0]])
        pts[0, 0] = 0.0
        c = SkyCoord(x=pts[0] * u.kpc, y=pts[1] * u.kpc, z=pts[2] * u.kpc,
                     v_x=np.zeros(4) * u.km / u.s, v_y=np.zeros(4) * u.km / u.s,
                     v_z=np.zeros(4) * u.km / u.s, frame=ICRS(), representation_type="cartesian",
                     differential_type="cartesian")
        g = c.transform_to(galcen_frame())
        xyz = np.array([g.x.to_value(u.kpc), g.y.to_value(u.kpc), g.z.to_value(u.kpc)])
        origin = xyz[:, 0]
        M = xyz[:, 1:] - origin[:, None]
        v_sun = np.array([g.v_x.to_value(u.km / u.s)[0], g.v_y.to_value(u.km / u.s)[0],
                          g.v_z.to_value(u.km / u.s)[0]])
        _affine.cache = (M, origin, v_sun)
    return _affine.cache


def to_galactocentric(ra, dec, d_kpc, pmra, pmdec, rv):
    """Arrays of position (kpc) and velocity (km/s), shape (3, n), in astropy's
    Galactocentric frame (x from the Sun's side toward the Galactic Centre,
    the Sun at x = -R0, rotation toward +y at the Sun). Uses the affine form
    of astropy's transformation (identical to ~1e-12; test_physics.py)."""
    M, origin, v_sun = _affine()
    a, dlt = np.radians(ra), np.radians(dec)
    d = np.asarray(d_kpc, float)
    ca, sa, cd, sd = np.cos(a), np.sin(a), np.cos(dlt), np.sin(dlt)
    rhat = np.array([cd * ca, cd * sa, sd * np.ones_like(ca)])
    ahat = np.array([-sa, ca, 0.0 * ca])
    dhat = np.array([-sd * ca, -sd * sa, cd * np.ones_like(ca)])
    rhat, ahat, dhat = (np.broadcast_to(h.reshape(3, -1), (3, d.size)) if h.ndim == 1 else h
                        for h in (rhat, ahat, dhat))
    x_icrs = rhat * d
    v_icrs = rhat * rv + K_KMS * d * (ahat * pmra + dhat * pmdec)
    return M @ x_icrs + origin[:, None], M @ v_icrs + v_sun[:, None]


def to_galactocentric_astropy(ra, dec, d_kpc, pmra, pmdec, rv):
    """Reference implementation through astropy (slow; used in tests)."""
    c = SkyCoord(ra=np.broadcast_to(ra, np.shape(d_kpc)) * u.deg,
                 dec=np.broadcast_to(dec, np.shape(d_kpc)) * u.deg,
                 distance=d_kpc * u.kpc, pm_ra_cosdec=pmra * u.mas / u.yr,
                 pm_dec=pmdec * u.mas / u.yr, radial_velocity=rv * u.km / u.s, frame=ICRS())
    g = c.transform_to(galcen_frame())
    x = np.array([g.x.to_value(u.kpc), g.y.to_value(u.kpc), g.z.to_value(u.kpc)])
    v = np.array([g.v_x.to_value(u.km / u.s), g.v_y.to_value(u.km / u.s),
                  g.v_z.to_value(u.km / u.s)])
    return x, v


def rotation_velocity(x, vcirc):
    """Circular disc velocity at positions x (3, n): speed vcirc (scalar,
    array or callable of R in kpc) along the direction of Galactic rotation."""
    R = np.hypot(x[0], x[1])
    vc = vcirc(R) if callable(vcirc) else vcirc
    return np.array([vc * x[1] / R, -vc * x[0] / R, np.zeros_like(R)])


KPC_PER_KMS_MYR = 1.0227121650537077e-3   # 1 km/s = 1.0227e-3 kpc/Myr


def launch(x, v, vcirc=V_CIRC_FLAT):
    """Straight-line flight back to the midplane (El-Badry et al. 2023, sec. 5).

    Returns t_disk (Myr; positive = crossed the plane in the past), the
    crossing point (3, n), and the ejection speed |v - v_rot(birth)|. When the
    crossing lies in the future (the star is still heading for the plane) the
    rotation velocity at the current position is subtracted, as they did."""
    t_disk = x[2] / (v[2] * KPC_PER_KMS_MYR)
    past = t_disk > 0
    xc = np.where(past, x - v * KPC_PER_KMS_MYR * t_disk, x)
    xc[2] = np.where(past, 0.0, x[2])
    v_ej = np.linalg.norm(v - rotation_velocity(xc, vcirc), axis=0)
    return t_disk, xc, v_ej


def posterior(row, rv, rv_err, ebv, rng, n_keep, prior="flat_MG", gauss=None,
              vcirc=V_CIRC_FLAT, v_max=V_EJECTION_MAX, zero_point=True):
    """Posterior samples for one star. ``row`` is its Gaia DR3 record."""
    plx_obs = float(row["parallax"])
    zp = parallax_zero_point(row) if zero_point else 0.0
    plx = plx_obs - zp
    s_plx = float(row["parallax_error"])
    s_a, s_d = float(row["pmra_error"]), float(row["pmdec_error"])
    r_pa, r_pd, r_ad = (float(row["parallax_pmra_corr"]), float(row["parallax_pmdec_corr"]),
                        float(row["pmra_pmdec_corr"]))
    cov = np.array([[s_plx ** 2, r_pa * s_plx * s_a, r_pd * s_plx * s_d],
                    [r_pa * s_plx * s_a, s_a ** 2, r_ad * s_a * s_d],
                    [r_pd * s_plx * s_d, r_ad * s_a * s_d, s_d ** 2]])
    # proper motions conditional on the parallax being 1/d
    gain = cov[1:, 0] / cov[0, 0]
    cov_cond = cov[1:, 1:] - np.outer(cov[1:, 0], cov[0, 1:]) / cov[0, 0]
    chol = np.linalg.cholesky(cov_cond)

    kept, tried, batch = [], 0, max(4 * n_keep, 20000)
    while sum(len(k["d"]) for k in kept) < n_keep and tried < 200 * n_keep:
        d = sample_distance(plx, s_plx, prior, batch, rng, gauss)
        resid = plx - 1.0 / d
        mu = (np.array([row["pmra"], row["pmdec"]], float)[:, None] - gain[:, None] * resid
              + chol @ rng.standard_normal((2, batch)))
        rvs = rv + rv_err * rng.standard_normal(batch)
        x, v = to_galactocentric(float(row["ra"]), float(row["dec"]), d, mu[0], mu[1], rvs)
        t_disk, xc, v_ej = launch(x, v, vcirc)
        ok = v_ej < v_max
        kept.append(dict(d=d[ok], pmra=mu[0][ok], pmdec=mu[1][ok], rv=rvs[ok], x=x[:, ok],
                         v=v[:, ok], t_disk=t_disk[ok], xc=xc[:, ok], v_ej=v_ej[ok]))
        tried += batch
    out = {k: (np.concatenate([b[k] for b in kept], axis=-1)) for k in kept[0]}
    n_ok = out["d"].size
    acceptance = n_ok / tried
    out = {k: (a[..., :n_keep]) for k, a in out.items()}
    d = out["d"]
    out["v_tan"] = K_KMS * np.hypot(out["pmra"], out["pmdec"]) * d
    out["v_tot"] = np.linalg.norm(out["v"], axis=0)
    out["v_helio"] = np.hypot(out["v_tan"], out["rv"])
    out["z"] = out["x"][2]
    ag = A_G_PER_EBV * ebv
    out["M_G0"] = float(row["phot_g_mean_mag"]) - ag - 5 * np.log10(100.0 * d)
    out["bp_rp0"] = float(row["bp_rp"]) - E_BPRP_PER_EBV * ebv if np.isfinite(row["bp_rp"]) else np.nan
    out["meta"] = dict(plx_corrected=plx, zero_point=zp, acceptance=acceptance, prior=prior,
                       ebv=ebv)
    return out


def summarize(a, q=(16, 50, 84)) -> dict:
    lo, med, hi = np.percentile(a, q)
    return dict(median=float(med), minus=float(med - lo), plus=float(hi - med),
                p2=float(np.percentile(a, 2.275)), p98=float(np.percentile(a, 97.725)))
