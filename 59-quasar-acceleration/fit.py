"""Weighted least-squares VSH fits to quasar proper motions, with Klioner et
al. (2021)'s outlier clipping, and three ways of putting error bars on them.

Each quasar contributes two equations (mu_ra*, mu_dec), weighted by the inverse
of its 2 x 2 covariance matrix from Gaia (the two errors and their correlation
rho). The fit is

    minimise  sum_k  r_k^T C_k^-1 r_k,   r_k = mu_k - B_k c

where B_k holds the VSH basis functions at the quasar's position (vsh.py) and c
the coefficients. The weighted residual of one quasar is X_k = sqrt(r^T C^-1 r)
(Klioner et al. 2021, Eq. 7); a quasar is an outlier when X_k exceeds kappa
times the median X of *all* quasars. Fit, flag, refit with the survivors, and
stop the first time a set of outliers repeats an earlier one (Lindegren 2018).

Error bars:

* ``formal``: the inverse normal matrix, as if Gaia's errors were exactly
  right, and ``scaled``: the same times the unit weight error squared (the
  reduced chi-square of the fit) -- Klioner et al.'s fit errors;
* ``sandwich``: N^-1 (sum_k g_k g_k^T) N^-1 with g_k = B_k^T C_k^-1 r_k, which
  takes the actual scatter of the residuals instead of the stated errors (the
  large-sample limit of resampling quasars with replacement);
* ``bootstrap``: resample the quasars with replacement and redo the whole fit,
  clipping included (Klioner et al.'s published errors are bootstrap errors).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

import numpy as np

import vsh
from config import A_G

CHUNK = 150_000


@dataclass
class Sample:
    """Positions (radians), proper motions and errors (uas/yr), correlation."""

    ra: np.ndarray
    dec: np.ndarray
    pmra: np.ndarray
    pmdec: np.ndarray
    pmra_error: np.ndarray
    pmdec_error: np.ndarray
    pmra_pmdec_corr: np.ndarray

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            setattr(self, name, np.ascontiguousarray(getattr(self, name), dtype=float))

    def __len__(self) -> int:
        return self.ra.size

    def subset(self, mask: np.ndarray) -> "Sample":
        return Sample(**{k: getattr(self, k)[mask] for k in self.__dataclass_fields__})

    def weights(self, use_corr: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Elements of the inverse covariance: w_aa, w_ad, w_dd."""
        rho = self.pmra_pmdec_corr if use_corr else np.zeros_like(self.ra)
        det = 1.0 - rho**2
        sa, sd = self.pmra_error, self.pmdec_error
        return 1.0 / (sa**2 * det), -rho / (sa * sd * det), 1.0 / (sd**2 * det)


@dataclass
class FitResult:
    lmax: int
    kappa: float
    term_list: list
    coef: np.ndarray
    cov_formal: np.ndarray
    uwe: float
    n_used: int
    outlier: np.ndarray
    iterations: int
    x: np.ndarray  # weighted residual X of every quasar
    cov_sandwich: np.ndarray | None = None
    cov_bootstrap: np.ndarray | None = None
    boot_coef: np.ndarray | None = None
    extra: dict = field(default_factory=dict)

    @property
    def cov_scaled(self) -> np.ndarray:
        return self.cov_formal * self.uwe**2

    def cov(self, kind: str = "scaled") -> np.ndarray:
        return {"formal": self.cov_formal, "scaled": self.cov_scaled,
                "sandwich": self.cov_sandwich, "bootstrap": self.cov_bootstrap}[kind]

    def vector(self, kind: str, frame: str = "icrs", cov: str = "scaled"):
        """Glide ('S') or rotation ('T') vector and its 3 x 3 covariance."""
        m = vsh.l1_vector_matrix(kind, self.term_list)
        if frame == "galactic":
            m = np.array(A_G) @ m
        c = self.cov(cov)
        return m @ self.coef, (None if c is None else m @ c @ m.T)

    @property
    def glide(self) -> np.ndarray:
        return self.vector("S")[0]

    @property
    def glide_cov(self) -> np.ndarray:
        return self.vector("S")[1]

    @property
    def rotation(self) -> np.ndarray:
        return self.vector("T")[0]

    @property
    def rotation_cov(self) -> np.ndarray:
        return self.vector("T")[1]

    @property
    def n_outliers(self) -> int:
        return int(self.outlier.sum())


def _chunks(n: int):
    for start in range(0, n, CHUNK):
        yield slice(start, min(start + CHUNK, n))


class Design:
    """The basis evaluated at every quasar, kept in memory as float32 when it
    fits (bootstraps reuse it hundreds of times), else recomputed per chunk."""

    def __init__(self, sample: Sample, term_list: list, max_bytes: float = 4e8) -> None:
        self.sample = sample
        self.term_list = term_list
        self.lmax = max(t.l for t in term_list)
        self.p = len(term_list)
        self.cached = 2 * len(sample) * self.p * 4 <= max_bytes
        if self.cached:
            self.b_ra = np.empty((len(sample), self.p), dtype=np.float32)
            self.b_dec = np.empty((len(sample), self.p), dtype=np.float32)
            for sl in _chunks(len(sample)):
                bra, bdec = vsh.basis(sample.ra[sl], sample.dec[sl], self.lmax, term_list)
                self.b_ra[sl], self.b_dec[sl] = bra, bdec

    def get(self, sl: slice) -> tuple[np.ndarray, np.ndarray]:
        if self.cached:
            return self.b_ra[sl].astype(float), self.b_dec[sl].astype(float)
        return vsh.basis(self.sample.ra[sl], self.sample.dec[sl], self.lmax, self.term_list)


def _normal_equations(design: Design, w, use: np.ndarray, mult: np.ndarray | None):
    s = design.sample
    w_aa, w_ad, w_dd = w
    p = design.p
    nmat = np.zeros((p, p))
    rhs = np.zeros(p)
    for sl in _chunks(len(s)):
        u = use[sl].astype(float)
        if mult is not None:
            u = u * mult[sl]
        keep = u > 0
        if not keep.any():
            continue
        bra, bdec = design.get(sl)
        bra, bdec, u = bra[keep], bdec[keep], u[keep]
        aa, ad, dd = (w_aa[sl][keep] * u, w_ad[sl][keep] * u, w_dd[sl][keep] * u)
        mra, mdec = s.pmra[sl][keep], s.pmdec[sl][keep]
        cross = bra.T @ (ad[:, None] * bdec)
        nmat += bra.T @ (aa[:, None] * bra) + cross + cross.T + bdec.T @ (dd[:, None] * bdec)
        rhs += bra.T @ (aa * mra + ad * mdec) + bdec.T @ (ad * mra + dd * mdec)
    return nmat, rhs


def _residual_x(design: Design, w, coef: np.ndarray) -> np.ndarray:
    s = design.sample
    w_aa, w_ad, w_dd = w
    x2 = np.empty(len(s))
    for sl in _chunks(len(s)):
        bra, bdec = design.get(sl)
        ra_r = s.pmra[sl] - bra @ coef
        de_r = s.pmdec[sl] - bdec @ coef
        x2[sl] = w_aa[sl] * ra_r**2 + 2 * w_ad[sl] * ra_r * de_r + w_dd[sl] * de_r**2
    return np.sqrt(np.maximum(x2, 0.0))


def _weighted_median(x: np.ndarray, mult: np.ndarray | None) -> float:
    if mult is None:
        return float(np.median(x))
    order = np.argsort(x)
    cw = np.cumsum(mult[order])
    return float(x[order][np.searchsorted(cw, 0.5 * cw[-1])])


def _sandwich(design: Design, w, coef: np.ndarray, use: np.ndarray, nmat_inv: np.ndarray):
    s = design.sample
    w_aa, w_ad, w_dd = w
    meat = np.zeros((design.p, design.p))
    for sl in _chunks(len(s)):
        keep = use[sl]
        if not keep.any():
            continue
        bra, bdec = design.get(sl)
        bra, bdec = bra[keep], bdec[keep]
        ra_r = s.pmra[sl][keep] - bra @ coef
        de_r = s.pmdec[sl][keep] - bdec @ coef
        a = w_aa[sl][keep] * ra_r + w_ad[sl][keep] * de_r
        d = w_ad[sl][keep] * ra_r + w_dd[sl][keep] * de_r
        g = bra * a[:, None] + bdec * d[:, None]
        meat += g.T @ g
    return nmat_inv @ meat @ nmat_inv


def _mask_key(mask: np.ndarray) -> str:
    return hashlib.sha1(np.packbits(mask).tobytes()).hexdigest()


def fit_design(design: Design, *, kappa: float = 3.0, use_corr: bool = True,
               mult: np.ndarray | None = None, max_iter: int = 50,
               sandwich: bool = True, weights=None) -> FitResult:
    """Fit with iterative clipping. ``mult`` are bootstrap multiplicities."""
    s = design.sample
    w = weights if weights is not None else s.weights(use_corr)
    present = np.ones(len(s), bool) if mult is None else mult > 0
    outlier = np.zeros(len(s), bool)
    seen = set()
    for it in range(1, max_iter + 1):
        nmat, rhs = _normal_equations(design, w, present & ~outlier, mult)
        nmat_inv = np.linalg.inv(nmat)
        coef = nmat_inv @ rhs
        x = _residual_x(design, w, coef)
        if not np.isfinite(kappa):
            new = np.zeros(len(s), bool)
        else:
            med = _weighted_median(x[present], None if mult is None else mult[present])
            new = present & (x > kappa * med)
        key = _mask_key(new)
        if key == _mask_key(outlier) or key in seen:
            if key != _mask_key(outlier):
                # the set cycles: refit with the repeated set, as Lindegren's rule implies
                outlier = new
                nmat, rhs = _normal_equations(design, w, present & ~outlier, mult)
                nmat_inv = np.linalg.inv(nmat)
                coef = nmat_inv @ rhs
                x = _residual_x(design, w, coef)
            break
        seen.add(_mask_key(outlier))
        outlier = new
    use = present & ~outlier
    if mult is None:
        n_used = int(use.sum())
        chi2 = float(np.sum(x[use] ** 2))
    else:
        n_used = int(mult[use].sum())
        chi2 = float(np.sum(mult[use] * x[use] ** 2))
    uwe = float(np.sqrt(chi2 / (2 * n_used - design.p)))
    res = FitResult(lmax=design.lmax, kappa=kappa, term_list=design.term_list, coef=coef,
                    cov_formal=nmat_inv, uwe=uwe, n_used=n_used, outlier=outlier,
                    iterations=it, x=x)
    if sandwich and mult is None:
        res.cov_sandwich = _sandwich(design, w, coef, use, nmat_inv)
    return res


def fit_vsh(sample: Sample, lmax: int = 3, kappa: float = 3.0, *, term_list=None,
            use_corr: bool = True, sandwich: bool = True, design: Design | None = None,
            weights=None) -> FitResult:
    tl = term_list if term_list is not None else vsh.terms(lmax)
    d = design if design is not None else Design(sample, tl)
    return fit_design(d, kappa=kappa, use_corr=use_corr, sandwich=sandwich, weights=weights)


def bootstrap(design: Design, n_rep: int, seed: int, *, kappa: float = 3.0,
              use_corr: bool = True, weights=None, progress: bool = False) -> np.ndarray:
    """Coefficients from ``n_rep`` resamplings of the quasars with replacement
    (multinomial counts), each refitted with clipping."""
    rng = np.random.default_rng(seed)
    n = len(design.sample)
    out = np.empty((n_rep, design.p))
    for r in range(n_rep):
        mult = rng.multinomial(n, np.full(n, 1.0 / n)).astype(float) if n < 50_000 else \
            np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        out[r] = fit_design(design, kappa=kappa, use_corr=use_corr, mult=mult,
                            sandwich=False, weights=weights).coef
        if progress and (r + 1) % 25 == 0:
            print(f"      bootstrap {r + 1}/{n_rep}", flush=True)
    return out
