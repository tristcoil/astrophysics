"""beta Pictoris b: an orbit seen edge-on, and its passages in front of the star.

    python betapic.py

Data: orbitize!'s example file (Nielsen et al. 2020, Table 1: relative
astrometry 2003-2018 from NaCo, NICI, GPI, SPHERE) and the CRIRES planet RV
of Snellen et al. (2014). The orbit is fitted with our own compiled
likelihood and nested sampling (fast), and cross-checked with orbitize!'s
MCMC (orbitize_fits.betapic).

Without the RV the sky motion cannot tell whether the planet passes in front
of the star in 2017 or behind it (the mirror orbit). The RV of -15.4 +- 1.7
km/s in Dec 2013 (planet approaching us while closing in on the star) picks
the branch: each posterior orbit and its mirror image are weighted by the RV
likelihood.

Outputs: conjunction epochs (front and behind) 1990-2070, the closest
projected approach and its distance in stellar radii and Hill radii.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

import config as C
import fastorbit
import kepler as K

PLX = (51.44, 0.12)          # mas, Hipparcos (van Leeuwen 2007), as in orbitize!'s beta Pic examples
MTOT = (1.75, 0.05)          # Msun, ditto
RV = (56643.0, -15.4, 1.7)   # MJD, planet RV relative to star (km/s), Snellen et al. 2014
#: angular diameter of beta Pic, 0.736 mas = 1.53 Rsun (Wang et al. 2016, Kervella et al. 2004)
STAR_DIAMETER_MAS = 0.736
#: planet mass for the Hill radius: Nielsen et al. 2020, 12.8 (+5.3 -3.2) MJup
MASS_B_MJUP = 12.8
LABELS = ["log10_a", "e", "cos_i", "varpi", "psi", "lambda", "plx", "mtot"]
PERIODIC = [3, 4, 5]


def load_data() -> pd.DataFrame:
    df = pd.read_csv(C.RAW / C.ORBITIZE_BETAPIC, comment="#")
    df = df[df.object == 1]
    ast = df[df.sep.notna()].copy()
    return ast.reset_index(drop=True)


class Model:
    def __init__(self):
        d = load_data()
        self.df = d
        self.epochs, self.eidx = np.unique(d.epoch.to_numpy(float), return_inverse=True)
        self.pidx = np.zeros(len(d), np.int64)
        self.is_seppa = np.ones(len(d), bool)
        self.obs = d[["sep", "pa"]].to_numpy(float)
        self.err = np.column_stack([d.sep_err, d.sep * np.radians(d.pa_err)]).astype(float)
        self.corr = np.zeros(len(d))
        self.jit = np.zeros(len(d))
        self.off = np.full(len(d), -1, np.int64)
        self.ndim = len(LABELS)

    def prior_transform(self, u):
        from scipy.special import ndtri
        x = np.empty(8)
        x[0] = 2.0 * u[0]                      # a log-uniform 1-100 au
        x[1] = 0.999 * u[1]
        x[2] = 1.0 - 2.0 * u[2]
        x[3:6] = 2 * np.pi * u[3:6]
        x[6] = PLX[0] + PLX[1] * ndtri(u[6])
        x[7] = MTOT[0] + MTOT[1] * ndtri(u[7])
        return x

    @staticmethod
    def elements(th):
        """-> (a, e, inc, aop, pan, tau) arrays, plx, mtot (vectorised)."""
        th = np.atleast_2d(th)
        a = 10 ** th[:, 0]
        e = th[:, 1]
        inc = np.arccos(np.clip(th[:, 2], -1, 1))
        aop = 0.5 * (th[:, 3] + th[:, 4])
        pan = np.mod(th[:, 3] - aop, 2 * np.pi)
        aop = np.mod(aop, 2 * np.pi)
        P = K.period_days(a, th[:, 7])
        M_ref = th[:, 5] - th[:, 3]
        tau = np.mod(np.mod((C.LAMBDA_REF_MJD - C.TAU_REF_MJD) / P, 1.0) - M_ref / (2 * np.pi), 1.0)
        return (a, e, inc, aop, pan, tau), th[:, 6], th[:, 7]

    def loglike(self, th):
        (a, e, inc, aop, pan, tau), plx, m = self.elements(th)
        if not (0 <= e[0] < 1):
            return -1e300
        el = np.array([[a[0], e[0], inc[0], aop[0], pan[0], tau[0]]])
        v = fastorbit.loglike(el, float(plx[0]), float(m[0]), np.zeros(1), C.TAU_REF_MJD, self.epochs,
                              self.pidx, self.eidx, self.is_seppa, self.obs, self.err, self.corr,
                              self.jit, self.off, np.zeros((0, 2)), np.zeros(0))
        return float(v) if np.isfinite(v) else -1e300


_M = None


def _ll(th):
    return _M.loglike(th)


def _pt(u):
    return _M.prior_transform(u)


def fit(nlive: int = 1000) -> np.ndarray:
    path = C.CACHE / "nested_betapic.npz"
    if path.exists():
        return np.load(path)["equal"]
    import dynesty
    global _M
    _M = Model()
    t0 = time.time()
    s = dynesty.NestedSampler(_ll, _pt, _M.ndim, nlive=nlive, bound="multi", sample="rslice",
                              periodic=PERIODIC, rstate=np.random.default_rng([C.RANDOM_SEED, 27321]))
    s.run_nested(dlogz=C.DLOGZ, print_progress=False)
    r = s.results
    eq = r.samples_equal(rstate=np.random.default_rng(C.RANDOM_SEED))
    np.savez_compressed(path, equal=eq, logz=r.logz[-1], logl=r.logl, samples=r.samples,
                        weights=r.importance_weights())
    print(f"beta Pic b: logZ = {r.logz[-1]:.2f}, {len(eq)} samples, {(time.time() - t0) / 60:.1f} min")
    return eq


def rv_branch(samples: np.ndarray, rng) -> tuple[tuple, np.ndarray, np.ndarray]:
    """Weight each orbit and its mirror image by the planet RV; draw one per sample.

    Returns the chosen elements, the planet RV (km/s) at the CRIRES epoch and
    the probability that the 'as sampled' branch is the right one.
    """
    (a, e, inc, aop, pan, tau), plx, m = Model.elements(samples)
    (_, _, _), (_, _, vz) = K.orbit_xyz(RV[0], a, e, inc, aop, pan, tau, m, velocity=True)
    kms = C.AU / C.YEAR / 1e3
    v = vz * kms                  # the mirror orbit has -v
    l1 = -0.5 * ((v - RV[1]) / RV[2]) ** 2
    l2 = -0.5 * ((-v - RV[1]) / RV[2]) ** 2
    p1 = 1.0 / (1.0 + np.exp(l2 - l1))
    keep = rng.uniform(size=len(p1)) < p1
    pan2 = np.where(keep, pan, np.mod(pan + np.pi, 2 * np.pi))
    aop2 = np.where(keep, aop, np.mod(aop + np.pi, 2 * np.pi))
    v2 = np.where(keep, v, -v)
    return (a, e, inc, aop2, pan2, tau), plx, m, v2, p1


def conjunctions(el, plx, m, t0_year=1990.0, t1_year=2070.0, dt_days=2.0):
    """Every local minimum of the projected separation, per posterior orbit.

    Returns a DataFrame: sample, year, mjd, sep_mas, sep_au, in_front (z < 0).
    """
    a, e, inc, aop, pan, tau = el
    t = Time_to_mjd(t0_year, t1_year, dt_days)
    rows = []
    for k in range(len(a)):
        x, y, z = K.orbit_xyz(t, a[k], e[k], inc[k], aop[k], pan[k], tau[k], m[k])
        s = np.hypot(x, y)
        i = np.where((s[1:-1] < s[:-2]) & (s[1:-1] < s[2:]))[0] + 1
        for j in i:
            # parabolic refinement of the minimum
            y0, y1, y2 = s[j - 1], s[j], s[j + 1]
            off = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            tm = t[j] + off * dt_days
            xm, ym, zm = K.orbit_xyz(tm, a[k], e[k], inc[k], aop[k], pan[k], tau[k], m[k])
            sm = float(np.hypot(xm, ym))
            rows.append(dict(sample=k, mjd=tm, sep_au=sm, sep_mas=sm * plx[k], in_front=bool(zm < 0)))
    out = pd.DataFrame(rows)
    from astropy.time import Time
    out["year"] = Time(out.mjd.to_numpy(), format="mjd").decimalyear
    return out


def Time_to_mjd(y0, y1, dt):
    from astropy.time import Time
    m0, m1 = Time([y0, y1], format="decimalyear").mjd
    return np.arange(m0, m1, dt)


def hill_radius_au(a, e, m_star, m_planet_mjup=MASS_B_MJUP):
    return a * (1 - e) * (m_planet_mjup * C.MJUP_MSUN / (3 * m_star)) ** (1 / 3)


def analyse(n: int = 3000) -> dict:
    eq = fit()
    rng = np.random.default_rng([C.RANDOM_SEED, 1])
    sub = eq[rng.choice(len(eq), min(n, len(eq)), replace=False)]
    el, plx, m, v_rv, p_branch = rv_branch(sub, rng)
    a, e, inc, aop, pan, tau = el
    P = K.period_days(a, m) / 365.25
    conj = conjunctions(el, plx, m)
    front17 = conj[(conj.year > 2016) & (conj.year < 2019) & conj.in_front]
    behind = conj[(conj.year > 2019) & (conj.year < 2040) & ~conj.in_front]
    nextfront = conj[(conj.year > 2030) & conj.in_front]
    nextfront = nextfront.sort_values("year").groupby("sample").head(1)
    rh = hill_radius_au(a, e, m)
    from astropy.time import Time

    def q(x):
        lo, med, hi = np.percentile(x, [16, 50, 84])
        return float(med), float(hi - med), float(med - lo)

    rstar_mas = STAR_DIAMETER_MAS / 2
    out = {
        "n_orbits": int(len(sub)),
        "a_au": q(a), "e": q(e), "i_deg": q(np.degrees(inc)), "P_yr": q(P),
        "Omega_deg": q(np.degrees(pan)), "omega_deg": q(np.degrees(aop)),
        "p_branch_front_2017": float(np.mean(p_branch)),
        "rv_predicted_kms": q(v_rv),
        "conj2017_mjd": q(front17.mjd), "conj2017_iso": Time(np.median(front17.mjd), format="mjd").iso[:10],
        "conj2017_sep_mas": q(front17.sep_mas), "conj2017_sep_au": q(front17.sep_au),
        "conj2017_sep_stellar_radii": q(front17.sep_mas / rstar_mas),
        "frac_planet_transits": float(np.mean(front17.sep_mas < rstar_mas)),
        "hill_radius_au": q(rh), "hill_radius_mas": q(rh * plx),
        "behind_year": q(behind.groupby("sample").year.first()),
        "next_front_year": q(nextfront.year),
        "next_front_sep_mas": q(nextfront.sep_mas),
    }
    # Hill sphere crossing in 2017: when the star is inside the planet's Hill radius
    t = Time_to_mjd(2016.5, 2018.8, 1.0)
    ins, outs = [], []
    for k in range(len(a)):
        x, y, z = K.orbit_xyz(t, a[k], e[k], inc[k], aop[k], pan[k], tau[k], m[k])
        inside = np.hypot(x, y) < rh[k]
        if inside.any() and np.mean(z[inside]) < 0:
            ins.append(t[inside][0])
            outs.append(t[inside][-1])
    out["hill_ingress_mjd"] = q(np.array(ins))
    out["hill_egress_mjd"] = q(np.array(outs))
    out["hill_ingress_iso"] = Time(np.median(ins), format="mjd").iso[:10]
    out["hill_egress_iso"] = Time(np.median(outs), format="mjd").iso[:10]
    out["_samples"] = dict(a=a, e=e, inc=inc, aop=aop, pan=pan, tau=tau, plx=plx, m=m)
    return out


if __name__ == "__main__":
    res = analyse()
    for k, v in res.items():
        if not k.startswith("_"):
            print(f"{k:28s} {v}")
