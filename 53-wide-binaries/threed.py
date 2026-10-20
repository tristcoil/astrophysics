"""3D test: 36 wide binaries with precise radial-velocity differences.

Chae et al. (2026, arXiv:2601.21728) compiled 36 pairs with both components'
radial velocities measured precisely enough (tens of m/s) to give the full 3D
relative velocity; they found a gravity boost gamma = G_eff/G ~ 1.6. Saad &
Ting (2026, arXiv:2603.11015) reanalysed the same data and got gamma = 1.00
(+0.24 -0.19) when the semi-major axis is a free parameter with a
population prior, but reproduced 1.56 when the 3D separation is instead a
geometric deprojection of the projected one.

We redo it independently and semi-analytically. For an isotropic population
of Kepler orbits with semi-major axis distribution P(a), eccentricity
distribution P(e) and uniform orbital phase, the density in phase space is
(Delaunay variables)

    f(x, v) = P(a) P(e) / (8 pi^3 mu^(3/2) a^(1/2) e),    mu = gamma G M,

with a and e computed from (x, v). Each pair's likelihood is f integrated
over its unknown line-of-sight separation z and over the measurement errors
of the three velocity components (and of the masses), by importance
sampling. Note that with P(a) ~ 1/a (Opik's law) the posterior over orbital
geometry given an exact projected separation is the geometric prior itself,
so "geometric deprojection" is this calculation with a log-flat P(a).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import AU_M, GM_SUN, RAW_DIR
import kinematics as KIN

KMS = 1e3
YR = 3.15576e7


def load():
    a = pd.read_csv(RAW_DIR / "chae2026_3d" / "chae_2026_data.csv")
    g = pd.read_csv(RAW_DIR / "chae2026_3d" / "chae_2026_gaia.csv")
    d = a.merge(g, on=["gaia_a", "gaia_b"], suffixes=("", "_g"), validate="one_to_one")
    return d


def observables(d, perspective: bool = True):
    """Projected separation vector (AU), sky relative velocity (km/s, 2D, B
    minus A, perspective removed with the systemic RV), its error, the RV
    difference and its error, masses."""
    ra1, de1, ra2, de2 = (d[c].to_numpy() for c in ("ra_a", "dec_a", "ra_b", "dec_b"))
    p1, p2 = d["parallax_a"].to_numpy(), d["parallax_b"].to_numpy()
    e1, e2 = d["parallax_error_a"].to_numpy(), d["parallax_error_b"].to_numpy()
    dist = KIN.weighted_distance(p1, e1, p2, e2)
    rv_sys = np.nanmean(np.stack([d["gaia_rv_a"].to_numpy(), d["gaia_rv_b"].to_numpy()]), axis=0)
    rv_sys = np.nan_to_num(rv_sys)
    vsky, rsky = KIN.sky_velocity(ra1, de1, d["pmra_a"].to_numpy(), d["pmdec_a"].to_numpy(), ra2, de2,
                                  d["pmra_b"].to_numpy(), d["pmdec_b"].to_numpy(), dist, dist, rv_sys,
                                  perspective=perspective)
    # star B minus star A, as 2D components in a sky basis at the midpoint
    n1, _, _ = KIN.unit(ra1, de1)
    n2, _, _ = KIN.unit(ra2, de2)
    nm = n1 + n2
    nm /= np.linalg.norm(nm, axis=1, keepdims=True)
    dr = (n2 - n1) * dist[:, None] * 206264.806            # AU, B minus A
    dr -= np.sum(dr * nm, axis=1, keepdims=True) * nm
    ex = dr / np.linalg.norm(dr, axis=1, keepdims=True)
    ey = np.cross(nm, ex)
    dv = -vsky                                              # sky_velocity gives A minus B
    v_par, v_perp = np.sum(dv * ex, axis=1), np.sum(dv * ey, axis=1)
    pm_err = 0.5 * np.sqrt(d["pmra_error_a"] ** 2 + d["pmdec_error_a"] ** 2 + d["pmra_error_b"] ** 2
                           + d["pmdec_error_b"] ** 2).to_numpy()
    v_err = KIN.K_PM * pm_err * dist / 1000.0
    # masses: Chae's (2026) values; 4 pairs have none, for which we use his
    # V-band mass-magnitude polynomial on the Gaia G magnitudes (Saad & Ting
    # set them to 1 Msun)
    import masses as ML
    m1, m2 = d["mass_a"].to_numpy(float), d["mass_b"].to_numpy(float)
    MG1 = d["phot_g_mean_mag_a"].to_numpy() - 5 * np.log10(dist / 10)
    MG2 = d["phot_g_mean_mag_b"].to_numpy() - 5 * np.log10(dist / 10)
    m1 = np.where(np.isfinite(m1), m1, ML.chae_mass(MG1))
    m2 = np.where(np.isfinite(m2), m2, ML.chae_mass(MG2))
    return dict(r_perp_au=np.linalg.norm(dr, axis=1), v_par=v_par, v_perp=v_perp, v_sky_err=v_err,
                v_z=d["vr_kms"].to_numpy(), v_z_err=d["vr_sigma_kms"].to_numpy(),
                m1=m1, m2=m2, dist=dist, chae_id=d["chae_id"].to_numpy())


def e_alpha_saad_ting(s_au):
    """Their separation-dependent eccentricity index (prior means)."""
    edges = np.array([0, 100, 300, 1000, 3000, 1e9])
    mu = np.array([-0.2, 0.4, 0.8, 1.0, 1.2])
    return mu[np.clip(np.digitize(s_au, edges) - 1, 0, 4)]


def log_pa(a_au, kind: str, a_break_au: float = 5000.0, a_min=10.0, a_max=1e6):
    """log P(a) (normalised on [a_min, a_max] AU)."""
    a = np.asarray(a_au, float)
    if kind == "opik":
        lp = -np.log(a) - np.log(np.log(a_max / a_min))
    elif kind == "bpl":     # a^-1 below the break, a^-1.6 above (Saad & Ting's baseline)
        grid = np.logspace(np.log10(a_min), np.log10(a_max), 4001)
        pdf = np.where(grid < a_break_au, 1 / grid, a_break_au**0.6 * grid**-1.6)
        norm = np.trapezoid(pdf, grid)
        lp = np.where(a < a_break_au, -np.log(a), 0.6 * np.log(a_break_au) - 1.6 * np.log(a)) - np.log(norm)
    else:
        raise ValueError(kind)
    return np.where((a > a_min) & (a < a_max), lp, -np.inf)


def pair_loglike(obs, i, gammas, *, pa: str = "bpl", n: int = 40000, seed: int = 0, rv_sys_kms: float = 0.04,
                 m_sigma: float = 0.05, vsys_sigma: float = 0.0, pm_jitter_mas: float = 0.0,
                 rv_student_nu: float | None = None, ignore_vr_angle: bool = False, vsys_in_rv: bool = False):
    """log L_i(gamma) on a grid, by importance sampling over z, the velocity
    errors and the masses."""
    rng = np.random.default_rng(seed + 1000 * i)
    rp = obs["r_perp_au"][i]
    # z: |z|/r_perp log-uniform over 1e-3 ... 1e3, random sign (proposal q)
    lt = rng.uniform(np.log(1e-3), np.log(1e3), n)
    z = rp * np.exp(lt) * np.where(rng.random(n) < 0.5, 1, -1)
    log_q_z = -np.log(np.abs(z)) - np.log(np.log(1e6)) - np.log(2)
    sv = obs["v_sky_err"][i]
    if vsys_sigma > 0:
        # an unconstrained systemic velocity (prior sigma vsys_sigma) enters only through perspective:
        # extra sky-velocity scatter of ~ vsys_sigma x angular separation (Saad & Ting's set-up)
        theta = rp / (obs["dist"][i] * 206264.806)
        sv = np.hypot(sv, vsys_sigma * theta)
    if pm_jitter_mas > 0:      # an extra proper-motion scatter, as in Saad & Ting's model
        sv = np.hypot(sv, KIN.K_PM * pm_jitter_mas * obs["dist"][i] / 1000.0)
    sz = np.hypot(obs["v_z_err"][i], rv_sys_kms)
    if vsys_sigma > 0 and vsys_in_rv:
        # the same unconstrained systemic velocity also shifts the RV difference, by its
        # component along the separation times the angular separation
        sz = np.hypot(sz, vsys_sigma * rp / (obs["dist"][i] * 206264.806))
    vx = obs["v_par"][i] + sv * rng.standard_normal(n)
    vy = obs["v_perp"][i] + sv * rng.standard_normal(n)
    zz = rng.standard_t(rv_student_nu, n) if rv_student_nu else rng.standard_normal(n)
    vz = obs["v_z"][i] + sz * zz
    if ignore_vr_angle:
        # use only the size of the projected separation, not its direction relative to the
        # sky velocity (as in a likelihood on |r_proj| alone): average over that angle
        ph = 2 * np.pi * rng.random(n)
        vx, vy = vx * np.cos(ph) - vy * np.sin(ph), vx * np.sin(ph) + vy * np.cos(ph)
    m = (obs["m1"][i] + obs["m2"][i]) * np.exp(m_sigma * rng.standard_normal(n))
    r = np.sqrt(rp**2 + z**2)                         # AU
    v2 = (vx**2 + vy**2 + vz**2) * (KMS * YR / AU_M) ** 2   # (AU/yr)^2
    # angular momentum |x cross v| with x = (rp, 0, z)
    hx, hy, hz = -z * vy, z * vx - rp * vz, rp * vy
    h2 = (hx**2 + hy**2 + hz**2) * (KMS * YR / AU_M) ** 2
    alpha_e = e_alpha_saad_ting(rp)
    out = np.empty(len(gammas))
    G_AU = 4 * np.pi**2 * GM_SUN / 1.32712440018e20     # AU^3 / (Msun yr^2) for GM_SUN in SI
    for k, gam in enumerate(gammas):
        mu = gam * G_AU * m
        inv_a = 2.0 / r - v2 / mu
        bound = inv_a > 0
        a = np.where(bound, 1.0 / np.where(bound, inv_a, 1.0), np.inf)
        e2 = np.clip(1.0 - h2 / (mu * a), 0.0, 1.0)
        e = np.sqrt(e2)
        lpe = np.log(alpha_e + 1.0) + alpha_e * np.log(np.maximum(e, 1e-12))
        lf = (log_pa(a, pa) + lpe - np.log(8 * np.pi**3) - 1.5 * np.log(mu) - 0.5 * np.log(a)
              - np.log(np.maximum(e, 1e-12)))
        lw = np.where(bound & (e < 1), lf - log_q_z, -np.inf)
        mx = np.max(lw)
        out[k] = mx + np.log(np.mean(np.exp(lw - mx))) if np.isfinite(mx) else -np.inf
    return out


def gamma_posterior(pa: str = "bpl", gammas=None, n: int = 40000, perspective: bool = True, drop=None,
                    vsys_sigma: float = 0.0, pm_jitter_mas: float = 0.0, rv_student_nu=None,
                    ignore_vr_angle: bool = False, vsys_in_rv: bool = False):
    """Posterior of gamma (flat in log gamma) from all pairs."""
    gammas = np.logspace(-0.6, 0.8, 141) if gammas is None else gammas
    obs = observables(load(), perspective=perspective)
    N = len(obs["r_perp_au"])
    keep = [i for i in range(N) if drop is None or i not in drop]
    ll = np.array([pair_loglike(obs, i, gammas, pa=pa, n=n, vsys_sigma=vsys_sigma, pm_jitter_mas=pm_jitter_mas,
                                rv_student_nu=rv_student_nu, ignore_vr_angle=ignore_vr_angle,
                                vsys_in_rv=vsys_in_rv) for i in keep])
    tot = ll.sum(0)
    post = np.exp(tot - tot.max())
    lg = np.log10(gammas)
    post /= np.trapezoid(post, lg)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(lg))])
    q = np.interp([0.16, 0.5, 0.84], cdf / cdf[-1], gammas)
    return dict(gammas=gammas, posterior_log10=post, q16=q[0], median=q[1], q84=q[2], per_pair=ll,
                n_pairs=len(keep))


if __name__ == "__main__":
    for pa in ("opik", "bpl"):
        r = gamma_posterior(pa)
        print(pa, f"gamma = {r['median']:.3f} (+{r['q84'] - r['median']:.3f} -{r['median'] - r['q16']:.3f})")


# --------------------------------------------------------------------------
# Allowing for pairs that are not bound Kepler orbits
# --------------------------------------------------------------------------

def spatial_density_perp(r_perp_au, pa: str, s_alpha, n_ae: int = 4000, n_z: int = 600, seed: int = 3):
    """n_perp(r_perp) = integral dz of the bound population's spatial density
    n(r) at r = sqrt(r_perp^2 + z^2), with the time-averaged radial density of
    a Kepler orbit p(r | a, e) = r / (pi a sqrt(a^2 e^2 - (r - a)^2))."""
    rng = np.random.default_rng(seed)
    # draw (a, e) from the priors
    grid = np.logspace(1, 6, 5001)
    pdf = np.exp(log_pa(grid, pa))
    cdf = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    a = np.interp(rng.random(n_ae), cdf / cdf[-1], grid)
    e = rng.random(n_ae) ** (1.0 / (s_alpha + 1.0))
    # z by Gauss-Legendre in u = asinh(z / r_perp) over a wide range
    x, w = np.polynomial.legendre.leggauss(n_z)
    umax = np.arcsinh(1e3)
    u = x * umax
    z = r_perp_au * np.sinh(u)
    dz = r_perp_au * np.cosh(u) * w * umax
    r = np.sqrt(r_perp_au**2 + z**2)
    R, A = np.meshgrid(r, a, indexing="ij")
    E = np.broadcast_to(e, R.shape)
    arg = A**2 * E**2 - (R - A) ** 2
    p = np.where(arg > 0, R / (np.pi * A * np.sqrt(np.maximum(arg, 1e-300))), 0.0)
    nr = p.mean(axis=1) / (4 * np.pi * r**2)
    return float(np.sum(nr * dz))


def gamma_posterior_mixture(pa: str = "bpl", gammas=None, n: int = 40000, v_max: float = 3.0,
                            f_grid=None, drop=None):
    """Each pair is a bound orbit with probability 1 - f, or (with
    probability f) an unbound or contaminated pair whose 3D relative velocity
    is uniform within v_max km/s; f has a flat prior and is marginalised."""
    gammas = np.logspace(-0.6, 0.8, 141) if gammas is None else gammas
    f_grid = np.linspace(0.0, 0.5, 51) if f_grid is None else f_grid
    obs = observables(load())
    N = len(obs["r_perp_au"])
    keep = [i for i in range(N) if drop is None or i not in drop]
    v_au_yr = KMS * YR / AU_M
    p_out = 1.0 / (4.0 / 3.0 * np.pi * (v_max * v_au_yr) ** 3)        # per (AU/yr)^3
    lb = []
    for i in keep:
        ll = pair_loglike(obs, i, gammas, pa=pa, n=n)
        nperp = spatial_density_perp(obs["r_perp_au"][i], pa, e_alpha_saad_ting(obs["r_perp_au"][i]))
        lb.append(ll - np.log(nperp))          # log p(v_obs | x_perp, gamma, bound)
    lb = np.array(lb)                           # (pairs, gammas)
    lo = np.log(p_out)
    # joint grid over (f, gamma)
    tot = np.zeros((len(f_grid), len(gammas)))
    for j, f in enumerate(f_grid):
        tot[j] = np.sum(np.logaddexp(np.log1p(-f) + lb, np.log(max(f, 1e-300)) + lo), axis=0)
    post2 = np.exp(tot - tot.max())
    lg = np.log10(gammas)
    post_g = post2.sum(0)
    post_g /= np.trapezoid(post_g, lg)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (post_g[1:] + post_g[:-1]) * np.diff(lg))])
    q = np.interp([0.16, 0.5, 0.84], cdf / cdf[-1], gammas)
    post_f = post2.sum(1)
    post_f /= post_f.sum()
    # per-pair probability of being an outlier at the posterior mode
    jm, km = np.unravel_index(np.argmax(tot), tot.shape)
    fo = f_grid[jm]
    w_out = np.exp(np.log(max(fo, 1e-300)) + lo - np.logaddexp(np.log1p(-fo) + lb[:, km],
                                                               np.log(max(fo, 1e-300)) + lo))
    return dict(gammas=gammas, posterior_log10=post_g, q16=q[0], median=q[1], q84=q[2], f_grid=f_grid,
                posterior_f=post_f, f_mode=float(fo), p_outlier=w_out, chae_id=obs["chae_id"][keep],
                n_pairs=len(keep), v_max=v_max)


def synthetic_obs(gamma_true: float, n_pairs: int, rng, *, pa: str = "bpl", r_range=(3000.0, 30000.0),
                  sv=0.012, sz=0.06, masses=(0.9, 0.7)):
    """A synthetic 36-pair-like sample: Kepler orbits from the population
    priors (P(a), the separation-dependent eccentricity law, isotropic
    orientation, uniform phase), kept if the projected separation lands in
    r_range, observed with Gaussian errors sv (sky) and sz (RV) in km/s."""
    grid = np.logspace(1, 6, 5001)
    pdf = np.exp(log_pa(grid, pa))
    cdf = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    out = {k: [] for k in ("r_perp_au", "v_par", "v_perp", "v_z", "m1", "m2")}
    G_AU = 4 * np.pi**2
    v_au = KMS * YR / AU_M
    while len(out["r_perp_au"]) < n_pairs:
        a = np.interp(rng.random(), cdf / cdf[-1], grid)
        al = e_alpha_saad_ting(a)
        e = rng.random() ** (1.0 / (al + 1.0))
        M = 2 * np.pi * rng.random()
        E = M
        for _ in range(60):
            E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        mu = gamma_true * G_AU * (masses[0] + masses[1])
        x = np.array([a * (np.cos(E) - e), a * np.sqrt(1 - e * e) * np.sin(E), 0.0])
        n_ = np.sqrt(mu / a**3)
        v = np.array([-a * n_ * np.sin(E), a * n_ * np.sqrt(1 - e * e) * np.cos(E), 0.0]) / (1 - e * np.cos(E))
        A = rng.standard_normal((3, 3))
        Q, R = np.linalg.qr(A)
        Q = Q * np.sign(np.diag(R))
        x, v = Q @ x, Q @ v
        rp = np.hypot(x[0], x[1])
        if not (r_range[0] < rp < r_range[1]):
            continue
        ex = x[:2] / rp
        ey = np.array([-ex[1], ex[0]])
        vk = v / v_au                                   # km/s
        out["r_perp_au"].append(rp)
        out["v_par"].append(vk[0] * ex[0] + vk[1] * ex[1] + sv * rng.standard_normal())
        out["v_perp"].append(vk[0] * ey[0] + vk[1] * ey[1] + sv * rng.standard_normal())
        out["v_z"].append(vk[2] + sz * rng.standard_normal())
        out["m1"].append(masses[0])
        out["m2"].append(masses[1])
    obs = {k: np.array(v) for k, v in out.items()}
    obs["v_sky_err"] = np.full(n_pairs, sv)
    obs["v_z_err"] = np.full(n_pairs, sz)
    obs["dist"] = np.full(n_pairs, 100.0)
    return obs


def injection_test(gamma_true: float = 1.0, n_sets: int = 8, n_pairs: int = 36, pa: str = "bpl", n: int = 20000):
    gammas = np.logspace(-0.6, 0.8, 141)
    meds = []
    for k in range(n_sets):
        rng = np.random.default_rng(900 + k)
        obs = synthetic_obs(gamma_true, n_pairs, rng, pa=pa)
        ll = np.array([pair_loglike(obs, i, gammas, pa=pa, n=n, rv_sys_kms=0.0, m_sigma=0.0)
                       for i in range(n_pairs)]).sum(0)
        post = np.exp(ll - ll.max())
        cdf = np.cumsum(post) / post.sum()
        meds.append(float(np.interp(0.5, cdf, gammas)))
    return np.array(meds)
