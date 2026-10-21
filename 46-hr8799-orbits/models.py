"""Orbit models: parameters, priors and likelihoods.

Four orbit models, defined as in Wang et al. (2018, Sec. 3), with their priors:

* ``unconstrained``      four independent orbits sharing parallax and star mass
* ``coplanar``           one shared inclination and node
* ``near1248``           period ratios uniform in 1.8-2.4 (b:c, c:d) and 1.8-2.2 (d:e)
* ``near1248_coplanar``  both, plus e < 0.2 for every planet

Priors (Wang et al.): a log-uniform in 1-100 au (planet e only in the ratio
models), e uniform, i ~ sin i, Omega, omega and tau uniform; parallax and star
mass Gaussian.

Parametrisation. Nine years of data cover a few to 15 per cent of these
orbits. Over such an arc the data pin down where each planet is on the sky and
how fast it moves, but hardly its distance along the line of sight or its
speed towards us. In orbital elements the allowed orbits then form a thin,
curved ridge (eccentricity, periastron and phase change together as the
unseen depth changes), and samplers lose their way on it: our first runs
converged on different narrow pieces of the ridge from one run to the next.
So we sample what the data measure. Each planet gets its sky position (x, y,
mas) and sky velocity (vx, vy, mas/yr) on MJD STATE_REF_MJD. Its line-of-sight
position and velocity (z, vz) are extra parameters in the non-coplanar models;
in the coplanar models they follow from the shared plane. The elements then
come from the state vector. The posterior becomes a compact blob, and Wang's
priors are carried over exactly with the phase-space Jacobian
(Delaunay variables):

    d^3x d^3v = (1/2) mu^(3/2) a^(1/2) e sin(i)  da de di dOmega domega dM
    d^2x d^2v = (1/2) mu e / sqrt(1 - e^2) cos(i)^2  da de dM domega   (sky
                projection of a planar orbit of inclination i)

so the prior density in the sampled coordinates is the elements' prior
divided by this Jacobian (times parallax factors for mas units). The sampled
coordinates have uniform priors over a box that contains every orbit
compatible with the data, and the evidence is exact: Z = V_box E_box[L pi].

Mirror orbits: without planet radial velocities, z -> -z gives the same sky
motion. The non-coplanar models sample z >= 0 only; the coplanar models
sample Omega in [0, 180 deg). Both halves have the same likelihood, so the
evidence is unchanged.

Two extensions for the 1998-2025 data set ("ALL"):

* planet masses (uniform 1-20 MJup) entering through the Jacobi geometry in
  kepler.system_offsets: the star wobbles around the inner planets' barycentre;
* per-instrument astrometric calibration: a north-angle rotation and a
  plate-scale error with Gaussian priors from each instrument's published
  calibration, marginalised analytically (they enter the model linearly), and
  a per-instrument jitter added in quadrature (log-uniform 0.01-20 mas).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.special import ndtri

import config as C
import fastorbit
import kepler as K

DEG = np.pi / 180.0
LOG_2PI = np.log(2.0 * np.pi)


# ---------------------------------------------------------------------------
# Data, prepared once for fast likelihood evaluations
# ---------------------------------------------------------------------------

@dataclass
class Data:
    df: pd.DataFrame
    epochs: np.ndarray            # unique MJDs
    eidx: np.ndarray              # row -> epoch index
    pidx: np.ndarray              # row -> planet index (C.INNER_TO_OUTER order)
    is_seppa: np.ndarray
    obs: np.ndarray               # (n, 2): (ra, dec) or (sep, pa)
    err: np.ndarray               # (n, 2) in mas: (ra, dec) or (sep, sep * pa_err_rad)
    corr: np.ndarray
    sources: list[str]
    sidx: np.ndarray              # row -> source index
    offset_sources: list[str] = field(default_factory=list)   # sources with calibration nuisances
    jitter_sources: list[str] = field(default_factory=list)

    @classmethod
    def from_frame(cls, df: pd.DataFrame, calibrate: bool = False) -> "Data":
        df = df.reset_index(drop=True)
        epochs, eidx = np.unique(df.mjd.to_numpy(), return_inverse=True)
        pidx = np.array([C.INNER_TO_OUTER.index(p) for p in df.planet])
        seppa = (df.kind == "seppa").to_numpy()
        obs = np.where(seppa[:, None], df[["sep", "pa"]].to_numpy(), df[["ra", "dec"]].to_numpy())
        err = np.where(seppa[:, None],
                       np.column_stack([df.sep_err, df.sep * np.radians(df.pa_err)]),
                       df[["ra_err", "dec_err"]].to_numpy())
        corr = np.where(seppa, 0.0, df["corr"].fillna(0.0).to_numpy())
        sources = sorted(df.source.unique())
        sidx = np.array([sources.index(s) for s in df.source])
        offs = [s for s in sources if calibrate and s in C.CALIBRATION]
        jit = [s for s in sources if calibrate and s not in ("NICMOS-S11", "GRAVITY-L19")]
        return cls(df, epochs, eidx, pidx, seppa, obs, err, corr, sources, sidx, offs, jit)


def sky_state_estimate(df: pd.DataFrame, planet: str, t_ref: float, window_yr: float = 8.0):
    """Model-free position and velocity of a planet on t_ref, from a weighted
    cubic fit to the data within +-window_yr. Returns (x, y, vx, vy) in mas,
    mas/yr and their standard errors (inflated by sqrt(chi2_nu) if > 1)."""
    d = df[(df.planet == planet) & (abs(df.mjd - t_ref) < window_yr * 365.25)]
    if len(d) < 6:
        d = df[df.planet == planet]
    t = (d.mjd.to_numpy() - t_ref) / 365.25
    A = np.column_stack([np.ones_like(t), t, t * t, t ** 3])
    out, err = [], []
    for col, ecol in (("ra", "ra_err"), ("dec", "dec_err")):
        w = 1.0 / d[ecol].to_numpy() ** 2
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        beta = cov @ (A.T @ (w * d[col].to_numpy()))
        res = d[col].to_numpy() - A @ beta
        chi2nu = float(np.sum(w * res ** 2) / max(len(t) - 4, 1))
        cov *= max(1.0, chi2nu)
        out.append(beta[:2])
        err.append(np.sqrt(np.diag(cov))[:2])
    (x, vx), (y, vy) = out
    (ex, evx), (ey, evy) = err
    return np.array([x, y, vx, vy]), np.array([ex, ey, evx, evy])


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

#: codes for fastorbit.loglike_state
KIND_CODE = {"unconstrained": 0, "coplanar": 1, "near1248": 2, "near1248_coplanar": 3}


class OrbitModel:
    """A parametrised orbit model with a prior transform and a likelihood."""

    def __init__(self, kind: str, data: Data, *, parallax=C.W18_PARALLAX, mstar=C.W18_MSTAR,
                 fit_masses: bool = False, mass_range=(1.0, 20.0), t_state=None):
        if kind not in C.MODELS:
            raise KeyError(kind)
        self.kind = kind
        self.data = data
        self.parallax = parallax
        self.mstar = mstar
        self.fit_masses = fit_masses
        self.mass_range = mass_range
        self.t_state = C.STATE_REF_MJD if t_state is None else t_state
        self.coplanar = kind in ("coplanar", "near1248_coplanar")
        self.ratios = kind in ("near1248", "near1248_coplanar")
        self.emax = C.LOW_E_MAX if kind == "near1248_coplanar" else C.ECC_RANGE[1]
        #: the sampled depths cover bound orbits with a < acap. The prior caps a at
        #: 100 au (planet e in the ratio models; the outer planets' a then reach ~4.8x
        #: a_e); every orbit the data allow is far below this.
        self.acap = 3.0 * C.SMA_RANGE[1]
        self._make_box()
        self.labels, self.periodic = self._layout()
        self.ndim = len(self.labels)
        self._prior_table()
        n_off = len(data.offset_sources)
        self._off_index = {s: data.sources.index(s) for s in data.offset_sources}
        if n_off:
            sig = np.array([[C.CALIBRATION[s][1], C.CALIBRATION[s][0] * DEG] for s in data.offset_sources])
            self._prior_prec = 1.0 / sig ** 2                    # (n_off, 2): scale, rotation
            self._logdet_prior = np.log(sig ** 2).sum(axis=1)
        else:
            self._prior_prec = np.zeros((0, 2))
            self._logdet_prior = np.zeros(0)
        self._jit_map = np.full(len(data.sources), -1)
        for k, s in enumerate(data.jitter_sources):
            self._jit_map[data.sources.index(s)] = k
        self._row_jit_idx = self._jit_map[data.sidx].astype(np.int64)
        q_of_source = np.full(len(data.sources), -1)
        for q, s in enumerate(data.offset_sources):
            q_of_source[data.sources.index(s)] = q
        self._off_row = q_of_source[data.sidx].astype(np.int64)
        self._idx = self._index_table()
        self._jit_idx = np.array([self.labels.index(f"log10_jit_{s}") for s in data.jitter_sources],
                                 dtype=np.int64)
        lo, hi = np.array(self._lo), np.array(self._hi)
        box = [k for k, n in enumerate(self.labels) if n.split("_")[0] in ("x", "y", "vx", "vy")
               or n in ("cos_i", "Omega")]
        self.log_vbox = float(np.sum(np.log(np.abs(hi[box] - lo[box]))))

    # -- the box of sky states ------------------------------------------------
    def _make_box(self):
        df = self.data.df
        plx_hi = self.parallax[0] + 6 * self.parallax[1]
        m_hi = self.mstar[0] + 6 * self.mstar[1] + (4 * self.mass_range[1] * C.MJUP_MSUN if self.fit_masses else 0)
        self.box = {}
        ests = {p: sky_state_estimate(df, p, self.t_state) for p in C.INNER_TO_OUTER}
        # With planet masses free, the sampled states are Jacobi coordinates: they
        # sit off the measured (star-relative) ones by the star's wobble around the
        # inner planets' barycentre, sum_j (m_j / M) X_j - up to ~30 mas for b at
        # 20 MJup each. Widen each box by the largest such shift.
        m_ratio = (self.mass_range[1] * C.MJUP_MSUN / max(self.mstar[0] - 3 * self.mstar[1], 0.5)
                   if self.fit_masses else 0.0)
        shift = np.zeros(4)
        for p in C.INNER_TO_OUTER:
            est, err = ests[p]
            # +-12 sigma of a model-free estimate, at least 15 mas and 2 mas/yr;
            # diagnostics.box_margins checks that no posterior comes near an edge
            half = np.maximum(12 * err, [15.0, 15.0, 2.0, 2.0]) + shift
            shift = shift + m_ratio * np.array([np.hypot(*est[:2])] * 2 + [np.hypot(*est[2:])] * 2)
            r_sky = max(np.hypot(est[0], est[1]) - half[0], 50.0) / plx_hi       # au, generous
            v_esc = np.sqrt(2 * C.G_AU_YR * m_hi / r_sky) * plx_hi                # mas/yr
            self.box[p] = dict(lo=est - half, hi=est + half,
                               z=(0.0, 2.0 * C.SMA_RANGE[1] * plx_hi), vz=(-1.05 * v_esc, 1.05 * v_esc))

    # -- parameter layout -------------------------------------------------
    def _layout(self):
        labels, periodic = [], []
        if self.coplanar:
            labels += ["cos_i", "Omega"]
            periodic.append(1)
        for p in C.INNER_TO_OUTER:                                      # e, d, c, b
            labels += [f"x_{p}", f"y_{p}", f"vx_{p}", f"vy_{p}"]
            if not self.coplanar:
                labels += [f"z_{p}", f"vz_{p}"]
        labels += ["plx", "mstar"]
        if self.fit_masses:
            labels += [f"m_{p}" for p in C.INNER_TO_OUTER]
        labels += [f"log10_jit_{s}" for s in self.data.jitter_sources]
        return labels, periodic

    # -- prior -------------------------------------------------------------
    def _prior_table(self):
        """Per parameter: (lo, hi) of a uniform prior, or (mu, sigma) of a Gaussian."""
        lo, hi, gauss = [], [], []
        for k, name in enumerate(self.labels):
            if name == "cos_i":
                a, b = 0.0, 1.0
            elif name == "Omega":
                a, b = 0.0, np.pi
            elif name.split("_")[0] in ("x", "y", "vx", "vy"):
                q, p = name.split("_")
                j = ("x", "y", "vx", "vy").index(q)
                a, b = self.box[p]["lo"][j], self.box[p]["hi"][j]
            elif name.startswith("z_"):
                a, b = self.box[name[2:]]["z"]
            elif name.startswith("vz_"):
                a, b = self.box[name[3:]]["vz"]
            elif name == "plx":
                a, b = self.parallax
                gauss.append(k)
            elif name == "mstar":
                a, b = self.mstar
                gauss.append(k)
            elif name.startswith("m_"):
                a, b = self.mass_range
            elif name.startswith("log10_jit_"):
                a, b = -2.0, np.log10(20.0)
            else:
                raise KeyError(name)
            lo.append(a)
            hi.append(b)
        self._lo, self._hi = np.array(lo, float), np.array(hi, float)
        self._gauss = np.array(gauss, int)
        self._unif = np.array([k for k in range(len(self.labels)) if k not in gauss], int)

    def prior_transform(self, u):
        x = self._lo + u * (self._hi - self._lo)
        g = self._gauss
        x[g] = self._lo[g] + self._hi[g] * ndtri(np.clip(u[g], 1e-15, 1 - 1e-15))
        if not self.coplanar:
            fastorbit.fill_depths(x, np.asarray(u, float), self._idx, self.fit_masses, C.MJUP_MSUN,
                                  C.G_AU_YR, self.acap)
        return x

    def log_prior_density(self, theta):
        """log density of the sampling prior (box x Gaussians), up to a constant."""
        th = np.asarray(theta, float)
        u = self._unif
        lo, hi = np.minimum(self._lo[u], self._hi[u]), np.maximum(self._lo[u], self._hi[u])
        per = np.isin(u, self.periodic)
        if not (per | ((th[u] >= lo) & (th[u] <= hi))).all():
            return -np.inf
        g = self._gauss
        return float(-0.5 * (((th[g] - self._lo[g]) / self._hi[g]) ** 2).sum())

    def _index_table(self):
        """Parameter positions for fastorbit.unpack_state (see its docstring)."""
        idx = np.full(40, -1, dtype=np.int64)
        pos = {n: k for k, n in enumerate(self.labels)}
        for k, p in enumerate(C.INNER_TO_OUTER):
            base = 8 * k
            for j, q in enumerate(("x", "y", "vx", "vy", "z", "vz")):
                idx[base + j] = pos.get(f"{q}_{p}", -1)
            if self.fit_masses:
                idx[base + 6] = pos[f"m_{p}"]
        if self.coplanar:
            idx[32], idx[33] = pos["cos_i"], pos["Omega"]
        idx[35], idx[36] = pos["plx"], pos["mstar"]
        return idx

    def _ratio_bounds(self):
        b = np.zeros((4, 2))
        for k, p in enumerate(C.INNER_TO_OUTER):
            if k > 0:
                b[k] = C.PERIOD_RATIO_RANGES[f"{p}/{C.INNER_TO_OUTER[k - 1]}"]
        return b

    # -- parameters -> physical ----------------------------------------------
    def physical(self, theta):
        """Elements (a, e, inc, aop, pan, tau) per planet, plx, M*, masses, jitters.

        Works on one parameter vector (shape ndim) or many (shape (n, ndim)).
        """
        th = np.asarray(theta, float)
        col = {n: th[..., k] for k, n in enumerate(self.labels)}
        plx, mstar = col["plx"], col["mstar"]
        masses = {p: (col[f"m_{p}"] * C.MJUP_MSUN if self.fit_masses else 0.0 * mstar)
                  for p in C.INNER_TO_OUTER}
        if self.coplanar:
            ci = col["cos_i"]
            ti = np.sqrt(np.maximum(1.0 - ci * ci, 0.0)) / ci
            cO, sO = np.cos(col["Omega"]), np.sin(col["Omega"])
        el, state = {}, {}
        m_inside = mstar
        for p in C.INNER_TO_OUTER:
            m_k = m_inside + masses[p]
            x, y, vx, vy = (col[f"{q}_{p}"] for q in ("x", "y", "vx", "vy"))
            if self.coplanar:
                z = ti * (cO * x - sO * y)
                vz = ti * (cO * vx - sO * vy)
            else:
                z, vz = col[f"z_{p}"], col[f"vz_{p}"]
            r = np.stack([x, y, z], axis=-1) / plx[..., None]
            v = np.stack([vx, vy, vz], axis=-1) / plx[..., None]
            a, e, inc, aop, pan, M = K.elements_from_state(r, v, C.G_AU_YR * m_k)
            P = K.period_days(np.where(a > 0, a, np.nan), m_k)
            tau = np.mod(np.mod((self.t_state - C.TAU_REF_MJD) / P, 1.0) - M / (2 * np.pi), 1.0)
            el[p] = (a, e, inc, aop, pan, tau)
            state[p] = (x, y, vx, vy, z, vz)
            m_inside = m_k
        jit = np.array([10.0 ** col[f"log10_jit_{s}"] for s in self.data.jitter_sources])
        return {"elements": el, "plx": plx, "mstar": mstar, "masses": masses, "jitter": jit,
                "state": state}

    # -- likelihood -------------------------------------------------------------
    def model_positions(self, phys, epochs=None):
        """(ra, dec) in mas for every planet at the data epochs: arrays (4, n_epochs)."""
        t = self.data.epochs if epochs is None else epochs
        sky = K.system_offsets(t, phys["elements"], phys["plx"], phys["mstar"],
                               phys["masses"] if self.fit_masses else None)
        ra = np.stack([sky[p][0] for p in C.INNER_TO_OUTER])
        dec = np.stack([sky[p][1] for p in C.INNER_TO_OUTER])
        return ra, dec

    def residuals(self, phys):
        """Residuals (data - model) in mas, and d(model)/d(scale, rotation)."""
        d = self.data
        ra, dec = self.model_positions(phys)
        mra, mdec = ra[d.pidx, d.eidx], dec[d.pidx, d.eidx]
        msep = np.hypot(mra, mdec)
        mpa = np.arctan2(mra, mdec)
        r = np.empty((len(d.pidx), 2))
        rd = ~d.is_seppa
        r[rd, 0] = d.obs[rd, 0] - mra[rd]
        r[rd, 1] = d.obs[rd, 1] - mdec[rd]
        sp = d.is_seppa
        r[sp, 0] = d.obs[sp, 0] - msep[sp]
        dpa = np.mod(np.radians(d.obs[sp, 1]) - mpa[sp] + np.pi, 2 * np.pi) - np.pi
        r[sp, 1] = dpa * d.obs[sp, 0]
        A = np.empty((len(d.pidx), 2, 2))
        A[rd, :, 0] = np.column_stack([mra[rd], mdec[rd]])
        A[rd, :, 1] = np.column_stack([mdec[rd], -mra[rd]])
        A[sp, 0, 0] = msep[sp]
        A[sp, 1, 0] = 0.0
        A[sp, 0, 1] = 0.0
        A[sp, 1, 1] = d.obs[sp, 0]
        return r, A, (mra, mdec)

    def loglike_phys(self, phys, return_offsets=False):
        """Data likelihood (numpy reference implementation of fastorbit.loglike)."""
        d = self.data
        r, A, _ = self.residuals(phys)
        s0, s1 = d.err[:, 0], d.err[:, 1]
        if len(d.jitter_sources):
            j = np.where(self._row_jit_idx >= 0, phys["jitter"][np.maximum(self._row_jit_idx, 0)], 0.0)
            s0 = np.sqrt(s0 ** 2 + j ** 2)
            s1 = np.sqrt(s1 ** 2 + j ** 2)
        c01 = d.corr * s0 * s1
        det = s0 ** 2 * s1 ** 2 - c01 ** 2
        i00, i11, i01 = s1 ** 2 / det, s0 ** 2 / det, -c01 / det
        chi2 = i00 * r[:, 0] ** 2 + i11 * r[:, 1] ** 2 + 2 * i01 * r[:, 0] * r[:, 1]
        lnl = -0.5 * (chi2.sum() + np.log(det).sum() + 2 * len(chi2) * LOG_2PI)
        offsets = {}
        for s, k in self._off_index.items():
            m = d.sidx == k
            Am, rm = A[m], r[m]
            NiA0 = i00[m, None] * Am[:, 0, :] + i01[m, None] * Am[:, 1, :]
            NiA1 = i01[m, None] * Am[:, 0, :] + i11[m, None] * Am[:, 1, :]
            F = (Am[:, 0, :, None] * NiA0[:, None, :] + Am[:, 1, :, None] * NiA1[:, None, :]).sum(axis=0)
            b = (NiA0 * rm[:, 0, None] + NiA1 * rm[:, 1, None]).sum(axis=0)
            q = list(d.offset_sources).index(s)
            Hm = np.diag(self._prior_prec[q]) + F
            sol = np.linalg.solve(Hm, b)
            lnl += 0.5 * b @ sol - 0.5 * (self._logdet_prior[q] + np.linalg.slogdet(Hm)[1])
            if return_offsets:
                offsets[s] = (sol, np.linalg.inv(Hm))
        return (lnl, offsets) if return_offsets else lnl

    def log_target_prior(self, phys) -> float:
        """log of Wang et al.'s prior density expressed in the sampled coordinates
        (numpy reference of the prior term inside fastorbit.loglike_state)."""
        el, plx = phys["elements"], float(phys["plx"])
        lp = 0.0
        if not self.coplanar:
            m_in = float(phys["mstar"])
            for p in C.INNER_TO_OUTER:
                m_in += float(phys["masses"][p])
                mu = C.G_AU_YR * m_in
                st = phys["state"][p]
                st = [float(np.asarray(v)) for v in st]
                zmax, rho, vs2 = fastorbit.bound_region(st[0], st[1], st[2], st[3], plx, mu, self.acap)
                vm = fastorbit.vz_max(st[4] / plx, rho, vs2, mu, self.acap) if zmax > 0 else 0.0
                if zmax <= 0 or vm <= 0:
                    return -np.inf
                lp += np.log(zmax * plx) + np.log(2 * vm * plx)
        m_inside = float(phys["mstar"])
        prev_P = None
        for k, p in enumerate(C.INNER_TO_OUTER):
            a, e, inc, aop, pan, tau = (float(np.asarray(v)) for v in el[p])
            m_k = m_inside + float(phys["masses"][p])
            mu = C.G_AU_YR * m_k
            P = float(K.period_days(a, m_k)) if a > 0 else np.nan
            if not (a > 0 and C.ECC_RANGE[0] <= e <= self.emax):
                return -np.inf
            if self.ratios and k > 0:
                lo, hi = C.PERIOD_RATIO_RANGES[f"{p}/{C.INNER_TO_OUTER[k - 1]}"]
                rr = P / prev_P
                if not lo <= rr <= hi:
                    return -np.inf
                lp += np.log(1.5 * rr / (a * (hi - lo)))
            else:
                if not C.SMA_RANGE[0] <= a <= C.SMA_RANGE[1]:
                    return -np.inf
                lp += -np.log(a * np.log(C.SMA_RANGE[1] / C.SMA_RANGE[0]))
            lp += -np.log(self.emax - C.ECC_RANGE[0])
            if self.coplanar:
                ci = np.cos(inc)
                lp += -2 * np.log(2 * np.pi) - (4 * np.log(plx) + 2 * np.log(abs(ci)) + np.log(mu)
                                                + np.log(e) - np.log(2.0) - 0.5 * np.log(1 - e * e))
            else:
                # the sin(i)/2 of the prior cancels the Jacobian's; log 2: z >= 0 mirror half
                lp += np.log(2.0) - 3 * np.log(2 * np.pi) - 6 * np.log(plx) - 1.5 * np.log(mu) \
                    - 0.5 * np.log(a) - np.log(e)
            prev_P = P
            m_inside = m_k
        if self.coplanar:
            lp += -np.log(2 * np.pi)          # shared (cos i, Omega): 1/(4 pi) x 2 mirrors
        return lp

    def loglike(self, theta):
        """ln L + ln(prior / sampling prior): what nested sampling integrates (compiled)."""
        d = self.data
        return fastorbit.loglike_state(
            np.asarray(theta, float), self._idx, KIND_CODE[self.kind], self.fit_masses,
            C.MJUP_MSUN, C.G_AU_YR, C.TAU_REF_MJD, self.t_state, self.emax, C.ECC_RANGE[0],
            C.SMA_RANGE[0], C.SMA_RANGE[1], self.acap, self._ratio_bounds(), self.log_vbox,
            self._jit_idx, self._row_jit_idx, d.epochs, d.pidx, d.eidx, d.is_seppa, d.obs, d.err, d.corr,
            self._off_row, self._prior_prec, self._logdet_prior)

    def loglike_reference(self, theta):
        """numpy version of loglike (for tests)."""
        phys = self.physical(theta)
        lt = self.log_target_prior(phys)
        if not np.isfinite(lt):
            return -1e300
        return float(self.loglike_phys(phys) + lt + self.log_vbox)

    def data_loglike(self, theta):
        """ln L of the data alone (no prior terms)."""
        return float(self.loglike_phys(self.physical(theta)))

    def chi2(self, theta):
        """Plain chi-square (no jitter, no offsets) - for Wang et al.'s goodness-of-fit numbers."""
        phys = self.physical(theta)
        r, _, _ = self.residuals(phys)
        d = self.data
        s0, s1 = d.err[:, 0], d.err[:, 1]
        c01 = d.corr * s0 * s1
        det = s0 ** 2 * s1 ** 2 - c01 ** 2
        return float(((s1 ** 2 * r[:, 0] ** 2 + s0 ** 2 * r[:, 1] ** 2 - 2 * c01 * r[:, 0] * r[:, 1]) / det).sum())


def derived(model: OrbitModel, samples: np.ndarray) -> pd.DataFrame:
    """Table of physical quantities for an array of posterior samples."""
    phys = model.physical(samples)
    out = {"plx": phys["plx"], "mstar": phys["mstar"]}
    for p in C.INNER_TO_OUTER:
        a, e, inc, aop, pan, tau = phys["elements"][p]
        out[f"a_{p}"] = a
        out[f"e_{p}"] = e
        out[f"i_{p}"] = np.degrees(inc)
        out[f"omega_{p}"] = np.degrees(aop)
        out[f"Omega_{p}"] = np.degrees(pan)
        out[f"tau_{p}"] = tau
        m_k = phys["mstar"] + sum(phys["masses"][q] for q in C.INNER_TO_OUTER[: C.INNER_TO_OUTER.index(p) + 1])
        out[f"P_{p}"] = K.period_days(a, m_k) / 365.25
        if model.fit_masses:
            out[f"m_{p}"] = phys["masses"][p] / C.MJUP_MSUN
    df = pd.DataFrame(out)
    for outer, inner in (("b", "c"), ("c", "d"), ("d", "e")):
        df[f"P_{outer}/P_{inner}"] = df[f"P_{outer}"] / df[f"P_{inner}"]
    return df
