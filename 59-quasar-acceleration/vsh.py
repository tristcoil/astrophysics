"""Vector spherical harmonics (VSH) on the sky, in the conventions of
Mignard & Klioner (2012, A&A 547, A59), which Klioner et al. (2021) used.

A proper-motion field on the sky -- an arrow at every point -- can be written as
a sum of two kinds of patterns, just as a scalar map can be written as a sum of
spherical harmonics:

* **toroidal** patterns T_lm swirl around points (no sources or sinks); the
  l = 1 toroidal terms are a rigid **rotation** of the whole sky.
* **spheroidal** patterns S_lm flow out of some points and into others; the
  l = 1 spheroidal terms are a **glide**: every arrow points along great
  circles from one point of the sky to the opposite one. An acceleration of
  the observer produces exactly a glide towards the acceleration.

Higher degrees l describe smaller-scale patterns (about 180/l degrees across);
in Gaia's quasar proper motions they soak up the instrument's large-scale
systematics, so they do not leak into the glide.

Definitions (Mignard & Klioner 2012, Eqs. 11-12, 30-31):

    Y_lm  = (-1)^m sqrt[(2l+1)/(4 pi) (l-m)!/(l+m)!] P_lm(sin dec) exp(i m ra)
    T_lm  = [dY/d(dec) e_ra - (1/cos dec) dY/d(ra) e_dec] / sqrt(l(l+1))
    S_lm  = [(1/cos dec) dY/d(ra) e_ra + dY/d(dec) e_dec] / sqrt(l(l+1))
    V     = sum_l [ t_l0 T_l0 + s_l0 S_l0
                    + 2 sum_{m>0} (t^R_lm T^R_lm - t^I_lm T^I_lm
                                   + s^R_lm S^R_lm - s^I_lm S^I_lm) ]

Here Y_lm is the standard (Condon-Shortley) spherical harmonic at colatitude
90 deg - dec and longitude ra, which is what scipy.special.sph_legendre_p
returns. Fits use the *real orthonormal* basis {T_l0, S_l0, sqrt(2) T^R_lm,
-sqrt(2) T^I_lm, sqrt(2) S^R_lm, -sqrt(2) S^I_lm}, whose coefficients map
one-to-one onto Mignard & Klioner's t and s (``to_mk``). The sign on the
imaginary parts makes each coefficient equal to +sqrt(2) times theirs.

The l = 1 terms are the rotation vector w and the glide vector g
(Klioner et al. 2021, Eq. 6; Mignard & Klioner 2012, Sect. 4.2):

    t_10 = sqrt(8 pi/3) w_z,  t^R_11 = -sqrt(4 pi/3) w_x,  t^I_11 = sqrt(4 pi/3) w_y
    s_10 = sqrt(8 pi/3) g_z,  s^R_11 = -sqrt(4 pi/3) g_x,  s^I_11 = sqrt(4 pi/3) g_y

so that (Eqs. 3-4 of Klioner et al. 2021)

    rotation: mu_ra* = -w_x sin(dec) cos(ra) - w_y sin(dec) sin(ra) + w_z cos(dec)
              mu_dec = +w_x sin(ra) - w_y cos(ra)
    glide:    mu_ra* = -g_x sin(ra) + g_y cos(ra)
              mu_dec = -g_x sin(dec) cos(ra) - g_y sin(dec) sin(ra) + g_z cos(dec)

test_physics.py checks the basis is orthonormal over the sphere and that these
l = 1 relations hold.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import sph_legendre_p


@dataclass(frozen=True)
class Term:
    """One real basis function: kind 'T' or 'S', degree l, order m, part."""

    kind: str  # "T" (toroidal) or "S" (spheroidal)
    l: int
    m: int
    part: str  # "0" for m = 0, "R" or "I" for m > 0

    @property
    def label(self) -> str:
        return f"{self.kind.lower()}{self.l}{self.m}" + ("" if self.part == "0" else self.part)

    @property
    def mk_factor(self) -> float:
        """Our coefficient = mk_factor x Mignard & Klioner's (t or s)."""
        return {"0": 1.0, "R": math.sqrt(2.0), "I": math.sqrt(2.0)}[self.part]


def terms(lmax: int, lmin: int = 1) -> list[Term]:
    """The real basis, ordered by degree, then order, toroidal before spheroidal.
    2 lmax (lmax + 2) functions for lmin = 1."""
    out = []
    for l in range(lmin, lmax + 1):
        for kind in ("T", "S"):
            out.append(Term(kind, l, 0, "0"))
        for m in range(1, l + 1):
            for kind in ("T", "S"):
                out.append(Term(kind, l, m, "R"))
                out.append(Term(kind, l, m, "I"))
    return out


def n_params(lmax: int) -> int:
    return 2 * lmax * (lmax + 2)


def basis(ra: np.ndarray, dec: np.ndarray, lmax: int,
          term_list: list[Term] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the real orthonormal VSH basis at positions (radians).

    Returns (B_ra, B_dec), each of shape (n_sources, n_terms): the e_ra and
    e_dec components of every basis function at every source.
    """
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    tl = terms(lmax) if term_list is None else term_list
    theta = 0.5 * np.pi - dec
    cosd = np.maximum(np.cos(dec), 1e-300)
    b_ra = np.empty((ra.size, len(tl)))
    b_dec = np.empty((ra.size, len(tl)))
    cache: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
    trig: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for j, t in enumerate(tl):
        key = (t.l, t.m)
        if key not in cache:
            p, dp_dtheta = sph_legendre_p(t.l, t.m, theta, diff_n=1)
            # dY/d(dec) = -dY/d(theta); (1/cos dec) dY/d(ra) = i m Y / cos dec
            cache[key] = (-dp_dtheta, t.m * p / cosd)
        d_ddec, m_over_cos = cache[key]
        norm = 1.0 / math.sqrt(t.l * (t.l + 1))
        if t.m == 0:
            if t.kind == "T":
                b_ra[:, j], b_dec[:, j] = norm * d_ddec, 0.0
            else:
                b_ra[:, j], b_dec[:, j] = 0.0, norm * d_ddec
            continue
        if t.m not in trig:
            trig[t.m] = (np.cos(t.m * ra), np.sin(t.m * ra))
        c, s = trig[t.m]
        f = math.sqrt(2.0) * norm * (1.0 if t.part == "R" else -1.0)
        # Real and imaginary parts of exp(i m ra) times the components above:
        #   S_lm: e_ra comp = i (m P/cos dec) e^{im ra}, e_dec comp = dY/d(dec) part
        #   T_lm: e_ra comp = dY/d(dec) part,  e_dec comp = -i (m P/cos dec) e^{im ra}
        if t.kind == "S":
            if t.part == "R":
                b_ra[:, j], b_dec[:, j] = -m_over_cos * s, d_ddec * c
            else:
                b_ra[:, j], b_dec[:, j] = m_over_cos * c, d_ddec * s
        else:
            if t.part == "R":
                b_ra[:, j], b_dec[:, j] = d_ddec * c, m_over_cos * s
            else:
                b_ra[:, j], b_dec[:, j] = d_ddec * s, -m_over_cos * c
        b_ra[:, j] *= f
        b_dec[:, j] *= f
    return b_ra, b_dec


# --------------------------------------------------------------------------
# l = 1: rotation and glide vectors
# --------------------------------------------------------------------------
_S8 = math.sqrt(8.0 * math.pi / 3.0)
_S4 = math.sqrt(4.0 * math.pi / 3.0)


def l1_vector_matrix(kind: str, term_list: list[Term]) -> np.ndarray:
    """Matrix M (3 x n_terms) with vector = M @ coefficients, for kind 'T'
    (rotation w) or 'S' (glide g), in the equatorial frame of the fit."""
    m = np.zeros((3, len(term_list)))
    for j, t in enumerate(term_list):
        if t.l != 1 or t.kind != kind:
            continue
        mk_per_ours = 1.0 / t.mk_factor  # Mignard & Klioner coefficient per ours
        if t.m == 0:
            m[2, j] = mk_per_ours / _S8          # z = s10 / sqrt(8 pi/3)
        elif t.part == "R":
            m[0, j] = -mk_per_ours / _S4         # x = -s11^R / sqrt(4 pi/3)
        else:
            m[1, j] = mk_per_ours / _S4          # y = +s11^I / sqrt(4 pi/3)
    return m


def glide_field(ra: np.ndarray, dec: np.ndarray, g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Proper motions (mu_ra*, mu_dec) produced by a glide vector g (Eq. 4)."""
    sa, ca, sd, cd = np.sin(ra), np.cos(ra), np.sin(dec), np.cos(dec)
    return (-g[0] * sa + g[1] * ca,
            -g[0] * sd * ca - g[1] * sd * sa + g[2] * cd)


def rotation_field(ra: np.ndarray, dec: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Proper motions (mu_ra*, mu_dec) produced by a rotation vector w."""
    sa, ca, sd, cd = np.sin(ra), np.cos(ra), np.sin(dec), np.cos(dec)
    return (-w[0] * sd * ca - w[1] * sd * sa + w[2] * cd,
            w[0] * sa - w[1] * ca)


def power_per_degree(coef: np.ndarray, term_list: list[Term]) -> dict[str, np.ndarray]:
    """Power of the fitted field in each degree, toroidal and spheroidal:
    P_l = sum of squared orthonormal coefficients = Mignard & Klioner's
    t_l0^2 + 2 sum_m (t^R^2 + t^I^2). Returned as the RMS proper motion over
    the sphere that each degree contributes, sqrt(P_l / 4 pi)."""
    lmax = max(t.l for t in term_list)
    out = {"T": np.zeros(lmax), "S": np.zeros(lmax)}
    for c, t in zip(coef, term_list):
        out[t.kind][t.l - 1] += c * c
    return {k: np.sqrt(v / (4 * math.pi)) for k, v in out.items()}
