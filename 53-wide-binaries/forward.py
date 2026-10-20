"""Forward model of the observed (r_sky, v~) distribution of wide binaries.

The model follows Banik et al. (2024, Section 3), rebuilt independently:

* **Wide binaries**: semi-major axes from a broken power law (slope -1.6
  below a_break, beta above), eccentricities P(e) = (gamma + 1) e^gamma,
  isotropic orientations, uniform orbital phase; total masses drawn from the
  observed sample; orbits from the alpha-libraries of orbits.py (alpha = 0
  Newton, 1 QUMOND with the Galactic field).
* **Undetected close companions (CBs)**: each star independently, with
  probability f_CB, is itself an unresolved pair with mass ratio q
  (P(q) ~ q^0.4, plus 4% exact twins), semi-major axis log-flat over 1.5 dex
  below k_CB times the wide orbit's a, thermal eccentricity. The hidden star
  adds mass the observer does not count (blended light raises the mass
  estimate by less than the true mass) and makes the photocentre wobble.
  Unlike Banik et al., who assume constant photocentre velocity above
  a_int = 3.2 AU and none below, we fit a straight line to the photocentre's
  track over Gaia DR3's 34-month baseline, so fast inner orbits average out
  as they do in the real proper motions.
* **Line-of-sight contamination**: a fraction f_LOS of unrelated pairs with
  dN/(dr_sky dv~) ~ v~ (Banik et al. eq. 20), or the empirical shape from
  our shifted (chance-alignment) catalogue.

Everything enters through importance weights on one big Monte Carlo sample:
for each alpha (library), each k_CB grid value and each of the four CB cases
(none, A, B, both) we store 4D histograms over (a bin, e bin, r_sky pixel,
v~ pixel). A likelihood call then only reweights the (a, e) bins and
interpolates linearly in alpha and log k_CB, which takes milliseconds.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numba import njit, prange

from config import A0_SI, AU_M, GM_SUN
import masses as ML
import orbits as OB

# Banik et al. (2024) table 2: pixels fixed in advance
R_EDGES_BANIK = np.concatenate([np.linspace(2, 5, 7)[:-1], np.linspace(5, 10, 6)[:-1],
                                np.linspace(10, 16, 4)[:-1], np.linspace(16, 22, 3)[:-1],
                                np.linspace(22, 30, 3)])
V_EDGES_BANIK = np.concatenate([np.linspace(0, 1.6, 21)[:-1], np.linspace(1.6, 2.6, 6)[:-1],
                                np.linspace(2.6, 5.0, 6)])
GAIA_BASELINE_YR = (34.0 / 12.0)          # DR3: 25 Jul 2014 - 28 May 2017
LOG_K_GRID = np.arange(-3.0, -0.2499, 0.25)   # log10 k_CB: 0.1% ... 56%
P_EQM = 0.04                               # twin fraction among CBs (Banik et al. app. A)
Q_INDEX = 0.4                              # P(q) ~ q^0.4 (Korntreff et al. 2012)
A_INT_DEX = 1.5
R_MIN_AU, R_MAX_AU = 50.0, 100_000.0       # orbits that crash/escape get zero weight
#: MOND radius of 1 Msun in AU
RM1_AU = np.sqrt(GM_SUN / A0_SI) / AU_M
GM_KMS2_AU = GM_SUN / AU_M / 1e6


@dataclass
class Binning:
    r_edges: np.ndarray = field(default_factory=lambda: R_EDGES_BANIK.copy())
    v_edges: np.ndarray = field(default_factory=lambda: V_EDGES_BANIK.copy())
    a_lo_kau: float = 1.0
    a_hi_kau: float = 100.0
    n_a: int = 30
    n_e: int = 25

    @property
    def a_edges(self):
        return np.logspace(np.log10(self.a_lo_kau), np.log10(self.a_hi_kau), self.n_a + 1)

    @property
    def e_edges(self):
        return np.linspace(0.0, 1.0, self.n_e + 1)


# --------------------------------------------------------------------------
# Close-binary photocentre motion over the Gaia baseline
# --------------------------------------------------------------------------

@njit(cache=True, fastmath=True)
def _kepler_E(M, e):
    E = M if e < 0.8 else np.pi
    for _ in range(40):
        f = E - e * np.sin(E) - M
        d = 1.0 - e * np.cos(E)
        dE = f / d
        E -= dE
        if abs(dE) < 1e-12:
            break
    return E


@njit(cache=True, fastmath=True)
def _photocentre_slope(a_au, mtot, e, R00, R01, R02, R10, R11, R12, M0, T_yr, n_t):
    """Least-squares slope (AU/yr, 2 sky components) of the relative orbit's
    projected position over [0, T]: the proper motion Gaia would measure for a
    unit photocentre lever arm. R** are the first two rows of the rotation
    from the orbital plane to the sky frame."""
    P = np.sqrt(a_au**3 / mtot)
    n = 2.0 * np.pi / P
    b = a_au * np.sqrt(1.0 - e * e)
    st = 0.0
    sx = 0.0
    sy = 0.0
    stt = 0.0
    stx = 0.0
    sty = 0.0
    for k in range(n_t):
        t = T_yr * (k + 0.5) / n_t
        E = _kepler_E((M0 + n * t) % (2.0 * np.pi), e)
        xo = a_au * (np.cos(E) - e)
        yo = b * np.sin(E)
        x = R00 * xo + R01 * yo
        y = R10 * xo + R11 * yo
        st += t
        sx += x
        sy += y
        stt += t * t
        stx += t * x
        sty += t * y
    den = n_t * stt - st * st
    return (n_t * stx - st * sx) / den, (n_t * sty - st * sy) / den


@njit(parallel=True, cache=True)
def cb_velocities(a_int_au, m_cb, e_int, rot, M0, lever, T_yr, n_t, out):
    """Photocentre proper-motion velocity (km/s, 2D) for each CB sample."""
    for s in prange(len(a_int_au)):
        if lever[s] == 0.0 or a_int_au[s] <= 0.0:
            out[s, 0] = 0.0
            out[s, 1] = 0.0
            continue
        vx, vy = _photocentre_slope(a_int_au[s], m_cb[s], e_int[s], rot[s, 0, 0], rot[s, 0, 1],
                                    rot[s, 0, 2], rot[s, 1, 0], rot[s, 1, 1], rot[s, 1, 2], M0[s],
                                    T_yr, n_t)
        out[s, 0] = lever[s] * vx * 4.740470463533348
        out[s, 1] = lever[s] * vy * 4.740470463533348


def random_rotations(rng, n):
    """Uniformly random 3x3 rotation matrices (QR of Gaussian matrices)."""
    A = rng.standard_normal((n, 3, 3))
    Q, Rr = np.linalg.qr(A)
    d = np.sign(np.diagonal(Rr, axis1=1, axis2=2))
    return Q * d[:, None, :]


# --------------------------------------------------------------------------
# Wide-binary states from an orbit library
# --------------------------------------------------------------------------

@njit(parallel=True, cache=True)
def wb_states(states, meta, log_a_grid, e_grid, n_cosi, a_kau, e, cosi, u_a, u_e, u_s, view,
              m_true, out):
    """For each sample: r~ = r_sky/a, v~ vector (2 sky components, relative
    to the true mass) and validity, from the library by stochastic
    interpolation in log(a/r_M) and e. ``view`` holds 3 unit vectors per
    sample (line of sight, sky e1, sky e2)."""
    na = len(log_a_grid)
    ne = len(e_grid)
    dla = log_a_grid[1] - log_a_grid[0]
    de = e_grid[1] - e_grid[0]
    ns = states.shape[3]
    for s in prange(len(a_kau)):
        rM = RM1 * np.sqrt(m_true[s])
        la = np.log10(a_kau[s] * 1000.0 / rM)
        fi = (la - log_a_grid[0]) / dla
        if fi < 0.0:
            fi = 0.0
        if fi > na - 1:
            fi = na - 1.0
        ia = int(fi)
        if ia < na - 1 and u_a[s] < fi - ia:
            ia += 1
        fe = (e[s] - e_grid[0]) / de
        if fe < 0.0:
            fe = 0.0
        if fe > ne - 1:
            fe = ne - 1.0
        ie = int(fe)
        if ie < ne - 1 and u_e[s] < fe - ie:
            ie += 1
        ic = int(cosi[s] * n_cosi)
        if ic >= n_cosi:
            ic = n_cosi - 1
        k = int(u_s[s] * ns)
        if k >= ns:
            k = ns - 1
        ag = 10.0 ** log_a_grid[ia]
        x0 = states[ia, ie, ic, k, 0]
        x1 = states[ia, ie, ic, k, 1]
        x2 = states[ia, ie, ic, k, 2]
        v0 = states[ia, ie, ic, k, 3]
        v1 = states[ia, ie, ic, k, 4]
        v2 = states[ia, ie, ic, k, 5]
        # projections on the sky basis (e1, e2)
        px = x0 * view[s, 1, 0] + x1 * view[s, 1, 1] + x2 * view[s, 1, 2]
        py = x0 * view[s, 2, 0] + x1 * view[s, 2, 1] + x2 * view[s, 2, 2]
        qx = v0 * view[s, 1, 0] + v1 * view[s, 1, 1] + v2 * view[s, 1, 2]
        qy = v0 * view[s, 2, 0] + v1 * view[s, 2, 1] + v2 * view[s, 2, 2]
        rs = np.sqrt(px * px + py * py)
        sq = np.sqrt(rs)
        out[s, 0] = rs / ag
        out[s, 1] = qx * sq
        out[s, 2] = qy * sq
        scale = a_kau[s] * 1000.0 / ag     # AU per r_M-unit for this sample
        ok = (meta[ia, ie, ic, 0] * scale > R_MIN) and (meta[ia, ie, ic, 1] * scale < R_MAX)
        out[s, 3] = 1.0 if ok else 0.0


RM1 = RM1_AU
R_MIN = R_MIN_AU
R_MAX = R_MAX_AU


def _hist4(ia, ie, ir, iv, w, H):
    """Add weights w into H[a, e, r, v] (samples with any index -1 skipped)."""
    na, ne, nr, nv = H.shape
    ok = (ia >= 0) & (ie >= 0) & (ir >= 0) & (iv >= 0) & (w > 0)
    flat = ((ia[ok] * ne + ie[ok]) * nr + ir[ok]) * nv + iv[ok]
    H += np.bincount(flat, weights=w[ok], minlength=H.size).reshape(H.shape).astype(H.dtype)


def _digitize(x, edges):
    i = np.searchsorted(edges, x, side="right") - 1
    i[(x < edges[0]) | (x >= edges[-1]) | ~np.isfinite(x)] = -1
    return i


# --------------------------------------------------------------------------
# Library of 4D histograms
# --------------------------------------------------------------------------

@dataclass
class ModelLibrary:
    binning: Binning
    alphas: np.ndarray
    log_k: np.ndarray
    H0: np.ndarray        # (n_alpha, na, ne, nr, nv)       no CB
    H1: np.ndarray        # (n_alpha, 3, nk, na, ne, nr, nv) CB on A, on B, on both
    N0: np.ndarray        # (n_alpha, na, ne) valid samples (all outcomes) per (a, e) bin
    N1: np.ndarray        # (n_alpha, 3, na, ne)
    n_samples: int
    meta: dict


def draw_population(mA, mB, n, binning: Binning, rng):
    """Common random numbers shared by every alpha, k_CB and CB case."""
    pick = rng.integers(0, len(mA), n)
    u = lambda: rng.random(n)
    pop = dict(mA=mA[pick], mB=mB[pick],
               a_kau=binning.a_lo_kau * (binning.a_hi_kau / binning.a_lo_kau) ** u(),
               e=u(), cosi=u(), u_a=u(), u_e=u(), u_s=u())
    los = rng.standard_normal((n, 3))
    los /= np.linalg.norm(los, axis=1, keepdims=True)
    t = rng.standard_normal((n, 3))
    e1 = t - np.sum(t * los, axis=1, keepdims=True) * los
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(los, e1)
    pop["view"] = np.ascontiguousarray(np.stack([los, e1, e2], axis=1))
    # close-binary draws for star A and B
    for w in ("A", "B"):
        twin = rng.random(n) < P_EQM
        q = np.where(twin, 1.0, rng.random(n) ** (1.0 / (1.0 + Q_INDEX)))
        pop[f"q{w}"] = q
        pop[f"w{w}"] = rng.random(n)                       # position in the 1.5-dex window
        pop[f"eint{w}"] = np.sqrt(rng.random(n))            # thermal: P(e) = 2e
        pop[f"rot{w}"] = np.ascontiguousarray(random_rotations(rng, n)[:, :2, :])
        pop[f"M0{w}"] = 2 * np.pi * rng.random(n)
    return pop


def cb_properties(m_est, q, table):
    """True mass, primary mass and |q~ - L~| lever for stars seen at m_est
    that are unresolved pairs of mass ratio q (blended light)."""
    lm = np.log10(np.clip(m_est, 10 ** table["logm"][0], 10 ** table["logm"][-1]))
    fi = (lm - table["logm"][0]) / (table["logm"][1] - table["logm"][0])
    fj = q * (len(table["q"]) - 1)
    i = np.clip(fi.astype(int), 0, len(table["logm"]) - 2)
    j = np.clip(fj.astype(int), 0, len(table["q"]) - 2)
    ti, tj = fi - i, fj - j

    def interp(arr):
        return ((1 - ti) * ((1 - tj) * arr[i, j] + tj * arr[i, j + 1])
                + ti * ((1 - tj) * arr[i + 1, j] + tj * arr[i + 1, j + 1]))

    mtrue = m_est * interp(table["mtrue_over_mest"])
    Lf = interp(table["L_frac"])
    qt = q / (1.0 + q)
    return mtrue, np.abs(qt - Lf)


def build_library(mA, mB, *, binning: Binning | None = None, n: int = 1_000_000,
                  alphas=None, log_k=None, seed: int = 20261009, quiet: bool = False,
                  mass_tag: str = "", cb_mode: str = "relative", stab_frac: float = 0.3,
                  vt_cap=None, noise_sampler=None) -> ModelLibrary:
    """All histograms for one observed mass distribution (mA, mB).

    ``cb_mode`` 'relative' (Banik et al.): inner semi-major axes log-flat over
    1.5 dex below k_CB times the wide orbit's a, so k_CB is relative.
    'absolute': log-flat over 1.5 dex below A_max = 10^log_k AU, independent
    of the wide orbit but capped at ``stab_frac`` times its a (dynamical
    stability); the 'log_k' grid then holds log10(A_max / AU).
    ``vt_cap(m_est, r_sky_kau, sig_vt)``: optional per-sample upper limit on the
    observed v~ (a selection such as El-Badry et al.'s 5 Msun cut).
    ``noise_sampler(r_sky_kau, rng)``: optional per-sample measurement error
    on each sky component of v~ (drawn from the data's own errors at that
    separation), added before the cap and the binning."""
    import time
    binning = binning or Binning()
    alphas = OB.ALPHA_GRID if alphas is None else np.asarray(alphas)
    log_k = LOG_K_GRID if log_k is None else np.asarray(log_k)
    rng = np.random.default_rng(seed)
    pop = draw_population(np.asarray(mA, float), np.asarray(mB, float), n, binning, rng)
    table = ML.blended_mass_table()
    ia = _digitize(pop["a_kau"], binning.a_edges)
    ie = _digitize(pop["e"], binning.e_edges)
    na, ne = binning.n_a, binning.n_e
    nr, nv = len(binning.r_edges) - 1, len(binning.v_edges) - 1
    nA, nK = len(alphas), len(log_k)
    H0 = np.zeros((nA, na, ne, nr, nv), np.float32)
    H1 = np.zeros((nA, 3, nK, na, ne, nr, nv), np.float32)
    N0 = np.zeros((nA, na, ne))
    N1 = np.zeros((nA, 3, na, ne))
    # masses per case
    mtA, levA = cb_properties(pop["mA"], pop["qA"], table)
    mtB, levB = cb_properties(pop["mB"], pop["qB"], table)
    m_est = pop["mA"] + pop["mB"]
    cases = {0: (pop["mA"], pop["mB"]), 1: (mtA, pop["mB"]), 2: (pop["mA"], mtB), 3: (mtA, mtB)}
    libs = {}
    t0 = time.time()
    zx, zy = rng.standard_normal(n), rng.standard_normal(n)   # unit noise, scaled per sample below
    # CB photocentre velocities depend on k_CB but not on alpha
    DV = np.zeros((nK, 2, n, 2), np.float32)
    for kk, lk in enumerate(log_k):
        k = 10.0 ** lk
        for jw, (wname, mt, lev) in enumerate((("A", mtA, levA), ("B", mtB, levB))):
            if cb_mode == "relative":
                a_int = k * 10.0 ** (-A_INT_DEX * pop[f"w{wname}"]) * pop["a_kau"] * 1000.0
            else:
                a_int = np.minimum(k * 10.0 ** (-A_INT_DEX * pop[f"w{wname}"]),
                                   stab_frac * pop["a_kau"] * 1000.0)
            outv = np.empty((n, 2))
            cb_velocities(a_int, mt, pop[f"eint{wname}"], pop[f"rot{wname}"], pop[f"M0{wname}"],
                          lev, GAIA_BASELINE_YR, 12, outv)
            DV[kk, jw] = outv
    if not quiet:
        print(f"    CB photocentre velocities for {nK} k values ({time.time() - t0:.0f} s)", flush=True)
    for ai, alpha in enumerate(alphas):
        lib = OB.build_library(float(alpha), quiet=True)
        st = np.ascontiguousarray(lib["states"])
        meta = np.ascontiguousarray(lib["meta"])
        for c, (ma, mb) in cases.items():
            m_true = ma + mb
            out = np.empty((n, 4))
            wb_states(st, meta, OB.LOG_A_GRID, OB.E_GRID, len(OB.COSI_GRID), pop["a_kau"], pop["e"],
                      pop["cosi"], pop["u_a"], pop["u_e"], pop["u_s"], pop["view"], m_true, out)
            valid = out[:, 3] > 0
            r_sky_kau = out[:, 0] * pop["a_kau"]
            vc_true = np.sqrt(GM_KMS2_AU * m_true / (r_sky_kau * 1000.0))
            vc_est = np.sqrt(GM_KMS2_AU * m_est / (r_sky_kau * 1000.0))
            vwb = out[:, 1:3] * vc_true[:, None]               # km/s
            sig_vt = (noise_sampler(r_sky_kau, np.random.default_rng(seed + 7)) if noise_sampler is not None
                      else np.zeros(n))
            if noise_sampler is not None:
                vwb = vwb + np.stack([zx * sig_vt * vc_est, zy * sig_vt * vc_est], axis=1)
            ir = _digitize(r_sky_kau, binning.r_edges)
            w = valid.astype(np.float64)
            if c == 0:
                np.add.at(N0[ai], (ia[valid], ie[valid]), 1.0)
                vt = np.linalg.norm(vwb, axis=1) / vc_est
                if vt_cap is not None:
                    vt = np.where(vt < vt_cap(m_est, r_sky_kau, sig_vt), vt, np.inf)
                iv = _digitize(vt, binning.v_edges)
                _hist4(ia, ie, ir, iv, w, H0[ai])
                continue
            np.add.at(N1[ai, c - 1], (ia[valid], ie[valid]), 1.0)
            libs[(ai, c)] = (vwb, vc_est, ir, w, sig_vt)
        if not quiet:
            print(f"    alpha {alpha:+.2f}: WB states done ({time.time() - t0:.0f} s)", flush=True)
        for kk in range(nK):
            for c in (1, 2, 3):
                vwb, vc_est, ir, w, sig_vt = libs[(ai, c)]
                vobs = vwb + (DV[kk, 0] if c in (1, 3) else 0.0) - (DV[kk, 1] if c in (2, 3) else 0.0)
                vt = np.linalg.norm(vobs, axis=1) / vc_est
                if vt_cap is not None:
                    vt = np.where(vt < vt_cap(m_est, r_sky_kau, sig_vt), vt, np.inf)
                iv = _digitize(vt, binning.v_edges)
                _hist4(ia, ie, ir, iv, w, H1[ai, c - 1, kk])
        libs.clear()
        if not quiet:
            print(f"    alpha {alpha:+.2f}: CB histograms done ({time.time() - t0:.0f} s)", flush=True)
    return ModelLibrary(binning=binning, alphas=np.asarray(alphas), log_k=np.asarray(log_k), H0=H0,
                        H1=H1, N0=N0, N1=N1, n_samples=n,
                        meta=dict(seed=seed, mass_tag=mass_tag, n_obs_masses=len(mA), cb_mode=cb_mode,
                                  stab_frac=stab_frac, vt_cap=vt_cap is not None,
                                  noise=noise_sampler is not None))


# --------------------------------------------------------------------------
# Likelihood
# --------------------------------------------------------------------------

def a_weights(binning: Binning, a_break, beta, slope_in=-1.6):
    """Bin masses of the broken power law P(a) ~ a^-1.6 (a < a_break),
    ~ a^beta beyond (continuous), over the bin masses of the log-uniform
    proposal: the importance weight of every sample in each a bin."""
    ed = binning.a_edges
    x0, x1 = ed[:-1], ed[1:]

    def F(lo, hi, p):
        return np.where(hi > lo, (hi ** (p + 1.0) - lo ** (p + 1.0)) / (p + 1.0), 0.0)

    inner = F(x0, np.minimum(x1, a_break), slope_in)
    outer = a_break ** (slope_in - beta) * F(np.maximum(x0, a_break), x1, beta)
    return (inner + outer) / np.log(x1 / x0)


def e_weights(binning: Binning, gamma):
    ed = binning.e_edges
    return np.diff(ed ** (gamma + 1.0)) / np.diff(ed)


def los_template(binning: Binning, kind: str = "banik", shifted=None):
    """Chance-alignment (r_sky, v~) shape, normalised to 1 over the pixels."""
    r, v = binning.r_edges, binning.v_edges
    if kind == "banik":
        T = np.outer(np.diff(r), np.diff(v**2) / 2.0)
    else:
        rs, vt = shifted
        T, _, _ = np.histogram2d(rs, vt, bins=[r, v])
        T = T + 1e-3
    return T / T.sum()


@dataclass
class Model:
    lib: ModelLibrary
    los: np.ndarray
    acceptance: np.ndarray | None = None    # per-pixel selection (e.g. R < 0.01), applied to bound pairs

    def bound_fraction_maps(self, theta):
        """Normalised bound-population pixel probabilities for parameters
        theta = (a_break, beta, gamma, f_CB, log10 k_CB, alpha)."""
        a_break, beta, gamma, f_cb, lk, alpha = theta
        L = self.lib
        W = np.outer(a_weights(L.binning, a_break, beta), e_weights(L.binning, gamma))
        # alpha interpolation
        A = L.alphas
        ja = int(np.clip(np.searchsorted(A, alpha) - 1, 0, len(A) - 2))
        ta = (alpha - A[ja]) / (A[ja + 1] - A[ja])
        K = L.log_k
        jk = int(np.clip(np.searchsorted(K, lk) - 1, 0, len(K) - 2))
        tk = (lk - K[jk]) / (K[jk + 1] - K[jk])
        fc = np.array([(1 - f_cb) ** 2, f_cb * (1 - f_cb), f_cb * (1 - f_cb), f_cb**2])
        out = 0.0
        for jj, wa in ((ja, 1 - ta), (ja + 1, ta)):
            p0 = np.tensordot(W, L.H0[jj], axes=([0, 1], [0, 1])) / np.sum(W * L.N0[jj])
            tot = fc[0] * p0
            for c in range(3):
                norm = np.sum(W * L.N1[jj, c])
                pk = ((1 - tk) * np.tensordot(W, L.H1[jj, c, jk], axes=([0, 1], [0, 1]))
                      + tk * np.tensordot(W, L.H1[jj, c, jk + 1], axes=([0, 1], [0, 1]))) / norm
                tot = tot + fc[c + 1] * pk
            out = out + wa * tot
        return out

    def pixel_probabilities(self, theta7):
        a_break, beta, gamma, f_cb, lk, f_los, alpha = theta7
        pb = self.bound_fraction_maps((a_break, beta, gamma, f_cb, lk, alpha))
        if self.acceptance is not None:
            pb = pb * self.acceptance
        s = pb.sum()
        pb = pb / s if s > 0 else pb
        return (1 - f_los) * pb + f_los * self.los


def binomial_lnlike(p, counts):
    """Banik et al.'s eq. (21-22): product of binomial pixel likelihoods."""
    from scipy.special import gammaln
    N = counts.sum()
    p = np.clip(p, 1e-300, 1 - 1e-15)
    lnC = gammaln(N + 1) - gammaln(counts + 1) - gammaln(N - counts + 1)
    return float(np.sum(lnC + counts * np.log(p) + (N - counts) * np.log1p(-p)))


# --------------------------------------------------------------------------
# Direct simulation (for mock catalogues / injection tests)
# --------------------------------------------------------------------------

@njit(parallel=True, cache=True)
def wb_states_vec(states, meta, log_a_grid, e_grid, n_cosi, a_kau, e, cosi, u_a, u_e, u_s, view,
                  m_true, out):
    """Like wb_states, plus the unit vector of the projected separation:
    out[:, 0] r~, [1:3] v~ vector (true mass), [3] valid, [4:6] separation
    direction in the sky basis."""
    na = len(log_a_grid)
    ne = len(e_grid)
    dla = log_a_grid[1] - log_a_grid[0]
    de = e_grid[1] - e_grid[0]
    ns = states.shape[3]
    for s in prange(len(a_kau)):
        rM = RM1 * np.sqrt(m_true[s])
        la = np.log10(a_kau[s] * 1000.0 / rM)
        fi = (la - log_a_grid[0]) / dla
        if fi < 0.0:
            fi = 0.0
        if fi > na - 1:
            fi = na - 1.0
        ia = int(fi)
        if ia < na - 1 and u_a[s] < fi - ia:
            ia += 1
        fe = (e[s] - e_grid[0]) / de
        if fe < 0.0:
            fe = 0.0
        if fe > ne - 1:
            fe = ne - 1.0
        ie = int(fe)
        if ie < ne - 1 and u_e[s] < fe - ie:
            ie += 1
        ic = int(cosi[s] * n_cosi)
        if ic >= n_cosi:
            ic = n_cosi - 1
        k = int(u_s[s] * ns)
        if k >= ns:
            k = ns - 1
        ag = 10.0 ** log_a_grid[ia]
        x0 = states[ia, ie, ic, k, 0]
        x1 = states[ia, ie, ic, k, 1]
        x2 = states[ia, ie, ic, k, 2]
        v0 = states[ia, ie, ic, k, 3]
        v1 = states[ia, ie, ic, k, 4]
        v2 = states[ia, ie, ic, k, 5]
        px = x0 * view[s, 1, 0] + x1 * view[s, 1, 1] + x2 * view[s, 1, 2]
        py = x0 * view[s, 2, 0] + x1 * view[s, 2, 1] + x2 * view[s, 2, 2]
        qx = v0 * view[s, 1, 0] + v1 * view[s, 1, 1] + v2 * view[s, 1, 2]
        qy = v0 * view[s, 2, 0] + v1 * view[s, 2, 1] + v2 * view[s, 2, 2]
        rs = np.sqrt(px * px + py * py)
        sq = np.sqrt(rs)
        out[s, 0] = rs / ag
        out[s, 1] = qx * sq
        out[s, 2] = qy * sq
        scale = a_kau[s] * 1000.0 / ag
        ok = (meta[ia, ie, ic, 0] * scale > R_MIN) and (meta[ia, ie, ic, 1] * scale < R_MAX)
        out[s, 3] = 1.0 if ok else 0.0
        out[s, 4] = px / rs
        out[s, 5] = py / rs


def sample_broken_power_law(rng, n, a_lo, a_hi, a_break, beta, slope_in=-1.6):
    """Inverse-CDF draws of a (kAU) from the broken power law on [a_lo, a_hi]."""
    grid = np.logspace(np.log10(a_lo), np.log10(a_hi), 4001)
    pdf = np.where(grid < a_break, grid ** slope_in, a_break ** (slope_in - beta) * grid ** beta)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, grid)


def simulate(theta7, mA, mB, n, rng, *, a_lo=0.1, a_hi=100.0, cb_mode="relative", stab_frac=0.3,
             alpha=None, e_alpha_fn=None):
    """Draw n model binaries (masses from the observed pairs mA, mB) at
    parameters theta7 = (a_break, beta, gamma, f_CB, log k_CB, f_LOS, alpha).
    Line-of-sight pairs are not simulated here (the mocks are of bound
    systems). Returns per-sample arrays: r_sky_kau, vt_vec (2D, in units of the
    Newtonian v_c at r_sky for the *estimated* masses), sep_dir (2D unit
    vector), m_est, m_true, valid, cases."""
    a_break, beta, gamma, f_cb, lk, f_los, al = theta7
    al = al if alpha is None else alpha
    pick = rng.integers(0, len(mA), n)
    ma, mb = np.asarray(mA, float)[pick], np.asarray(mB, float)[pick]
    a_kau = sample_broken_power_law(rng, n, a_lo, a_hi, a_break, beta)
    if e_alpha_fn is None:
        e = rng.random(n) ** (1.0 / (gamma + 1.0))
    else:   # separation-dependent power law, e.g. Hwang et al.'s alpha(log10 a/AU)
        e = rng.random(n) ** (1.0 / (e_alpha_fn(a_kau * 1000.0) + 1.0))
    u = lambda: rng.random(n)
    los = rng.standard_normal((n, 3))
    los /= np.linalg.norm(los, axis=1, keepdims=True)
    t = rng.standard_normal((n, 3))
    e1 = t - np.sum(t * los, axis=1, keepdims=True) * los
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    view = np.ascontiguousarray(np.stack([los, e1, np.cross(los, e1)], axis=1))
    table = ML.blended_mass_table()
    hasA, hasB = rng.random(n) < f_cb, rng.random(n) < f_cb
    mt, dv = {}, {}
    for w, m, has in (("A", ma, hasA), ("B", mb, hasB)):
        twin = rng.random(n) < P_EQM
        q = np.where(twin, 1.0, rng.random(n) ** (1.0 / (1.0 + Q_INDEX)))
        mtrue, lev = cb_properties(m, q, table)
        mt[w] = np.where(has, mtrue, m)
        k = 10.0 ** lk
        wpos = rng.random(n)
        if cb_mode == "relative":
            a_int = k * 10.0 ** (-A_INT_DEX * wpos) * a_kau * 1000.0
        else:
            a_int = np.minimum(k * 10.0 ** (-A_INT_DEX * wpos), stab_frac * a_kau * 1000.0)
        outv = np.empty((n, 2))
        cb_velocities(a_int, mt[w], np.sqrt(rng.random(n)),
                      np.ascontiguousarray(random_rotations(rng, n)[:, :2, :]), 2 * np.pi * rng.random(n),
                      np.where(has, lev, 0.0), GAIA_BASELINE_YR, 12, outv)
        dv[w] = outv
    m_true = mt["A"] + mt["B"]
    m_est = ma + mb
    # alpha may fall between library values: pick one of the two neighbours
    A = OB.ALPHA_GRID
    j = int(np.clip(np.searchsorted(A, al) - 1, 0, len(A) - 2))
    tsel = (al - A[j]) / (A[j + 1] - A[j])
    use_hi = rng.random(n) < tsel
    out = np.empty((n, 6))
    for jj, sel in ((j, ~use_hi), (j + 1, use_hi)):
        if not np.any(sel):
            continue
        lib = OB.build_library(float(A[jj]), quiet=True)
        o = np.empty((int(sel.sum()), 6))
        wb_states_vec(np.ascontiguousarray(lib["states"]), np.ascontiguousarray(lib["meta"]), OB.LOG_A_GRID,
                      OB.E_GRID, len(OB.COSI_GRID), a_kau[sel], e[sel], u()[sel], u()[sel], u()[sel], u()[sel],
                      view[sel], m_true[sel], o)
        out[sel] = o
    r_sky = out[:, 0] * a_kau
    vc_true = np.sqrt(GM_KMS2_AU * m_true / (r_sky * 1000.0))
    vc_est = np.sqrt(GM_KMS2_AU * m_est / (r_sky * 1000.0))
    v = out[:, 1:3] * vc_true[:, None] + dv["A"] - dv["B"]
    return dict(r_sky_kau=r_sky, vt_vec=v / vc_est[:, None], sep_dir=out[:, 4:6], m_est=m_est,
                m_true=m_true, valid=out[:, 3] > 0, pick=pick, hasA=hasA, hasB=hasB, a_kau=a_kau, e=e)
