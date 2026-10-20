"""The Standard Model's couplings, run from the top-quark mass to the Planck scale.

How strongly particles interact depends on the energy at which you look:
quantum fluctuations screen or anti-screen each charge, and the couplings
"run" with the energy scale mu according to the renormalization group
equations (RGEs), d(coupling)/d ln mu = beta(couplings). The Higgs field's
self-coupling lambda is dragged down by the heavy top quark and, with the
measured masses, turns slightly negative at ~10^10 GeV: the Higgs potential
then dips below the value it has in our vacuum, which is therefore not the
lowest state available -- it is metastable.

Two steps, both from Buttazzo, Degrassi, Giardino, Giudice, Sala, Salvio &
Strumia (2013), "Investigating the near-criticality of the Higgs boson",
JHEP 12, 089 (arXiv:1307.3536v4, whose TeX source we read the formulas from):

1. **Matching** at mu = M_t: the MS-bar couplings lambda, y_t, g2, gY follow
   from the measured pole masses (M_h, M_t, M_W, M_Z) and the Fermi constant.
   Tree level, the full analytic one-loop corrections (their appendix A), and
   their numerical two-loop fits (section 2) are implemented; y_t also gets the
   three-loop pure-QCD term, g3 comes from alpha_s(M_Z) (their section 3.5).
2. **Running** with the beta functions at one, two or three loops (their
   appendix B, which collects Mihaila-Salomon-Steinhauser, Chetyrkin-Zoller
   and Bednyakov et al.), plus the four-loop QCD term for g3.

The Higgs effective potential at large field values, V = lambda_eff(h) h^4/4,
is implemented at one and two loops (their appendix C), and the vacuum decay
rate follows their section 6.3.

The state vector holds squared couplings and is integrated in s = ln mu^2,
exactly the form in which appendix B prints the equations:
    X = (g1^2, g2^2, g3^2, yt^2, yb^2, ytau^2, lambda, Gamma)
with g1 = sqrt(5/3) gY in SU(5) normalisation and Gamma = int gamma d ln mu
the integrated Higgs-field anomalous dimension. The potential is
V = -m^2|H|^2/2 + lambda |H|^4, so that M_h^2 = 2 lambda v^2 at tree level.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import quad
from scipy.special import spence

K = 1.0 / (16.0 * np.pi**2)
ZETA3 = 1.2020569031595942
NAMES = ("g1sq", "g2sq", "g3sq", "ytsq", "ybsq", "ytausq", "lam", "Gamma")

# Buttazzo et al.'s reference point for their fit formulas (their table 1)
REF = {"Mh": 125.15, "Mt": 173.34, "MW": 80.384, "MW_err": 0.014, "MZ": 91.1876,
       "alpha_s": 0.1184, "V": 246.21971}


# ---------------------------------------------------------------------------
# Beta functions (Buttazzo et al. 2013, appendix B), d X / d ln mu^2
# ---------------------------------------------------------------------------


def betas(X: np.ndarray, loops: int = 3) -> np.ndarray:
    """dX/ds for the state X (shape (8,) or (8, N)); loops = 1, 2 or 3.

    At three loops the pure-QCD four-loop term of g3 is included, as in the
    paper's numerical code. Bottom and tau Yukawas run at two loops.
    """
    G1, G2, G3, YT, YB, YTA, L = X[:7]
    k, k2, k3 = K, K * K, K**3

    # ---- one loop
    dG1 = k * G1**2 * (41 / 10)
    dG2 = k * G2**2 * (-19 / 6)
    dG3 = k * G3**2 * (-7.0)
    dYT = k * YT * (9 / 2 * YT + 3 / 2 * YB + YTA - 8 * G3 - 9 / 4 * G2 - 17 / 20 * G1)
    dYB = k * YB * (3 / 2 * YT + 9 / 2 * YB + YTA - 8 * G3 - 9 / 4 * G2 - 1 / 4 * G1)
    dYTA = k * YTA * (3 * YT + 3 * YB + 5 / 2 * YTA - 9 / 4 * G2 - 9 / 4 * G1)
    dL = k * (L * (12 * L + 6 * YT + 6 * YB + 2 * YTA - 9 / 2 * G2 - 9 / 10 * G1)
              - 3 * YT**2 - 3 * YB**2 - YTA**2
              + 9 / 16 * G2**2 + 27 / 400 * G1**2 + 9 / 40 * G2 * G1)
    gamma = k * (9 / 4 * G2 + 9 / 20 * G1 - 3 * YT - 3 * YB - YTA)

    if loops >= 2:
        dG1 = dG1 + k2 * G1**2 * (44 / 5 * G3 + 27 / 10 * G2 + 199 / 50 * G1
                                  - 17 / 10 * YT - 1 / 2 * YB - 3 / 2 * YTA)
        dG2 = dG2 + k2 * G2**2 * (12 * G3 + 35 / 6 * G2 + 9 / 10 * G1
                                  - 3 / 2 * YT - 3 / 2 * YB - 1 / 2 * YTA)
        dG3 = dG3 + k2 * G3**2 * (-26 * G3 + 9 / 2 * G2 + 11 / 10 * G1 - 2 * YT - 2 * YB)
        dYT = dYT + k2 * YT * (
            YT * (-12 * YT - 11 / 4 * YB - 9 / 4 * YTA - 12 * L + 36 * G3
                  + 225 / 16 * G2 + 393 / 80 * G1)
            + YB * (-1 / 4 * YB + 5 / 4 * YTA + 4 * G3 + 99 / 16 * G2 + 7 / 80 * G1)
            + YTA * (-9 / 4 * YTA + 15 / 8 * G2 + 15 / 8 * G1)
            + 6 * L**2 - 108 * G3**2 - 23 / 4 * G2**2 + 1187 / 600 * G1**2
            + 9 * G3 * G2 + 19 / 15 * G3 * G1 - 9 / 20 * G2 * G1)
        dYB = dYB + k2 * YB * (
            YT * (-1 / 4 * YT - 11 / 4 * YB + 5 / 4 * YTA + 4 * G3 + 99 / 16 * G2
                  + 91 / 80 * G1)
            + YB * (-12 * YB - 9 / 4 * YTA - 12 * L + 36 * G3 + 225 / 16 * G2
                    + 237 / 80 * G1)
            + YTA * (-9 / 4 * YTA + 15 / 8 * G2 + 15 / 8 * G1)
            + 6 * L**2 - 108 * G3**2 - 23 / 4 * G2**2 - 127 / 600 * G1**2
            + 9 * G3 * G2 + 31 / 15 * G3 * G1 - 27 / 20 * G2 * G1)
        dYTA = dYTA + k2 * YTA * (
            6 * L**2 - 23 / 4 * G2**2 + 1371 / 200 * G1**2 + 27 / 20 * G2 * G1
            + YT * (-27 / 4 * YT + 3 / 2 * YB - 27 / 4 * YTA + 20 * G3 + 45 / 8 * G2
                    + 17 / 8 * G1)
            + YB * (-27 / 4 * YB - 27 / 4 * YTA + 20 * G3 + 45 / 8 * G2 + 5 / 8 * G1)
            + YTA * (-3 * YTA - 12 * L + 165 / 16 * G2 + 537 / 80 * G1))
        dL = dL + k2 * (
            L**2 * (-156 * L - 72 * YT - 72 * YB - 24 * YTA + 54 * G2 + 54 / 5 * G1)
            + L * YT * (-3 / 2 * YT - 21 * YB + 40 * G3 + 45 / 4 * G2 + 17 / 4 * G1)
            + L * YB * (-3 / 2 * YB + 40 * G3 + 45 / 4 * G2 + 5 / 4 * G1)
            + L * YTA * (-1 / 2 * YTA + 15 / 4 * G2 + 15 / 4 * G1)
            + L * (-73 / 16 * G2**2 + 1887 / 400 * G1**2 + 117 / 40 * G2 * G1)
            + YT**2 * (15 * YT - 3 * YB - 16 * G3 - 4 / 5 * G1)
            + YT * (-9 / 8 * G2**2 - 171 / 200 * G1**2 + 63 / 20 * G2 * G1)
            + YB**2 * (-3 * YT + 15 * YB - 16 * G3 + 2 / 5 * G1)
            + YB * (-9 / 8 * G2**2 + 9 / 40 * G1**2 + 27 / 20 * G2 * G1)
            + YTA**2 * (5 * YTA - 6 / 5 * G1)
            + YTA * (-3 / 8 * G2**2 - 9 / 8 * G1**2 + 33 / 20 * G2 * G1)
            + 305 / 32 * G2**3 - 3411 / 4000 * G1**3 - 289 / 160 * G2**2 * G1
            - 1677 / 800 * G2 * G1**2)
        gamma = gamma + k2 * (
            YT * (-3 / 2 * YB - 17 / 8 * G1 - 45 / 8 * G2 - 20 * G3 + 27 / 4 * YT)
            - YTA * (15 / 8 * G1 + 15 / 8 * G2 - 9 / 4 * YTA)
            + YB * (-5 / 8 * G1 - 45 / 8 * G2 - 20 * G3 + 27 / 4 * YB)
            - 1293 / 800 * G1**2 - 27 / 80 * G2 * G1 + 271 / 32 * G2**2 - 6 * L**2)

    if loops >= 3:
        dG1 = dG1 + k3 * G1**2 * (
            YT * (189 / 16 * YT - 29 / 5 * G3 - 471 / 32 * G2 - 2827 / 800 * G1)
            + L * (-9 / 5 * L + 9 / 10 * G2 + 27 / 50 * G1)
            + 297 / 5 * G3**2 + 789 / 64 * G2**2 - 388613 / 24000 * G1**2
            - 3 / 5 * G3 * G2 - 137 / 75 * G3 * G1 + 123 / 160 * G2 * G1)
        dG2 = dG2 + k3 * G2**2 * (
            YT * (147 / 16 * YT - 7 * G3 - 729 / 32 * G2 - 593 / 160 * G1)
            + L * (-3 * L + 3 / 2 * G2 + 3 / 10 * G1)
            + 81 * G3**2 + 324953 / 1728 * G2**2 - 5597 / 1600 * G1**2
            + 39 * G3 * G2 - 1 / 5 * G3 * G1 + 873 / 160 * G2 * G1)
        dG3 = dG3 + k3 * G3**2 * (
            YT * (15 * YT - 40 * G3 - 93 / 8 * G2 - 101 / 40 * G1)
            + 65 / 2 * G3**2 + 109 / 8 * G2**2 - 523 / 120 * G1**2
            + 21 * G3 * G2 + 77 / 15 * G3 * G1 - 3 / 40 * G2 * G1) \
            + K**4 * G3**5 * (-2472.28)
        dYT = dYT + k3 * YT * (
            YT**2 * (58.6028 * YT + 198 * L - 157 * G3 - 1593 / 16 * G2 - 2437 / 80 * G1)
            + L * YT * (15 / 4 * L + 16 * G3 - 135 / 2 * G2 - 127 / 10 * G1)
            + YT * (363.764 * G3**2 + 16.990 * G2**2 - 24.422 * G1**2 + 48.370 * G3 * G2
                    + 18.074 * G3 * G1 + 34.829 * G2 * G1)
            + L**2 * (-36 * L + 45 * G2 + 9 * G1)
            + L * (-171 / 16 * G2**2 - 1089 / 400 * G1**2 + 117 / 40 * G2 * G1)
            - 619.35 * G3**3 + 169.829 * G2**3 + 16.099 * G1**3 + 73.654 * G3**2 * G2
            - 15.096 * G3**2 * G1 - 21.072 * G3 * G2**2 - 22.319 * G3 * G1**2
            - 321 / 20 * G3 * G2 * G1 - 4.743 * G2**2 * G1 - 4.442 * G2 * G1**2)
        dL = dL + k3 * (
            L**3 * (6011.35 * L + 873 * YT - 387.452 * G2 - 77.490 * G1)
            + L**2 * YT * (1768.26 * YT + 160.77 * G3 - 359.539 * G2 - 63.869 * G1)
            + L**2 * (-790.28 * G2**2 - 185.532 * G1**2 - 316.64 * G2 * G1)
            + L * YT**2 * (-223.382 * YT - 662.866 * G3 - 5.470 * G2 - 21.015 * G1)
            + L * YT * (356.968 * G3**2 - 319.664 * G2**2 - 74.8599 * G1**2
                        + 15.1443 * G3 * G2 + 17.454 * G3 * G1 + 5.615 * G2 * G1)
            + L * G2**2 * (-57.144 * G3 + 865.483 * G2 + 79.638 * G1)
            + L * G1**2 * (-8.381 * G3 + 61.753 * G2 + 28.168 * G1)
            + YT**3 * (-243.149 * YT + 250.494 * G3 + 74.138 * G2 + 33.930 * G1)
            + YT**2 * (-50.201 * G3**2 + 15.884 * G2**2 + 15.948 * G1**2 + 13.349 * G3 * G2
                       + 17.570 * G3 * G1 - 70.356 * G2 * G1)
            + YT * G3 * (16.464 * G2**2 + 1.016 * G1**2 + 11.386 * G2 * G1)
            + YT * G2**2 * (62.500 * G2 + 13.041 * G1)
            + YT * G1**2 * (10.627 * G2 + 11.117 * G1)
            + G3 * (7.536 * G2**3 + 0.663 * G1**3 + 1.507 * G2**2 * G1 + 1.105 * G2 * G1**2)
            - 114.091 * G2**4 - 1.508 * G1**4 - 37.889 * G2**3 * G1 + 6.500 * G2**2 * G1**2
            - 1.543 * G2 * G1**3)
        gamma = gamma + k3 * (
            -9 * G1 * L**2 - 45 * G2 * L**2 + 1.07 * G1**2 * L + 3.57 * G2 * G1 * L
            + 8.92 * G2**2 * L + 14.99 * G1**2 * YT + 14.13 * G1 * YT**2
            - 13.21 * G2 * G1 * YT - 8.73 * G3 * G1 * YT + 40.11 * G2 * YT**2
            + 79.05 * G3 * YT**2 + 23.40 * G2**2 * YT - 178.48 * G3**2 * YT
            - 7.57 * G2 * G3 * YT - 5.26 * G1**3 + 1.93 * G2 * G1**2 + 4.19 * G3 * G1**2
            + 1.81 * G2**2 * G1 - 158.51 * G2**3 + 28.57 * G2**2 * G3 + 36 * L**3
            + 135 / 2 * L**2 * YT - 45.00 * L * YT**2 - 60.13 * YT**3)

    # Gamma = int gamma d ln mu  ->  dGamma/d ln mu^2 = gamma / 2
    return np.array([dG1, dG2, dG3, dYT, dYB, dYTA, dL, 0.5 * gamma])


def run(X0: np.ndarray, mu0: float, mu1: float, *, loops: int = 3,
        step: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the RGEs from mu0 to mu1 (GeV) with fourth-order Runge-Kutta.

    X0 may be (8,) or (8, N) for N trajectories at once. Returns (mu, X) at
    every step, X with shape (n_steps + 1, 8[, N]). The step is in ln mu^2;
    0.05 keeps the integration error below 1e-10 (checked in test_physics.py
    against an adaptive DOP853 solution).
    """
    s0, s1 = 2 * np.log(mu0), 2 * np.log(mu1)
    n = max(1, int(np.ceil(abs(s1 - s0) / step)))
    h = (s1 - s0) / n
    X = np.array(X0, dtype=float)
    out = np.empty((n + 1, *X.shape))
    out[0] = X
    for i in range(n):
        k1 = betas(X, loops)
        k2 = betas(X + 0.5 * h * k1, loops)
        k3 = betas(X + 0.5 * h * k2, loops)
        k4 = betas(X + h * k3, loops)
        X = X + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = X
    mu = np.exp(0.5 * (s0 + h * np.arange(n + 1)))
    return mu, out


# ---------------------------------------------------------------------------
# Matching at mu = M_t (Buttazzo et al. 2013, sections 2-3 and appendix A)
# ---------------------------------------------------------------------------


def A0(M: float, mu: float) -> float:
    """Finite part of the one-point Passarino-Veltman function."""
    return 0.0 if M == 0 else M**2 * (1 - np.log(M**2 / mu**2))


def _xlogx(z):
    return 0.0 if z == 0 else z * np.log(z)


def B0(p: float, M1: float, M2: float, mu: float) -> float:
    """Real part of the finite two-point function,
    B0 = -int_0^1 ln |x M1^2 + (1-x) M2^2 - x(1-x) p^2| / mu^2 dx.

    The argument is p^2 (x - r1)(x - r2), and int_0^1 ln(x - r) dx =
    (1-r) ln(1-r) + r ln(-r) - 1 for real or complex r (real parts taken), so
    B0 has a closed form that stays finite through thresholds, p = M1 + M2,
    where the two roots merge."""
    if M1 == 0 or M2 == 0:            # one root sits at an end point
        M = M1 + M2
        return (2 - M**2 / p**2 * np.log(M**2 / mu**2)
                + (M**2 - p**2) / p**2 * np.log(abs(M**2 - p**2) / mu**2))
    a, b, c = p**2, M1**2 - M2**2 - p**2, M2**2      # f(x) = a x^2 + b x + c
    if a < 1e-10 * max(M1, M2) ** 2:                  # p -> 0: no roots in [0, 1]
        f = lambda x: -np.log((x * M1**2 + (1 - x) * M2**2 - x * (1 - x) * a) / mu**2)
        return quad(f, 0.0, 1.0, epsabs=1e-13, epsrel=1e-12)[0]
    sq = np.sqrt(complex(b * b - 4 * a * c))
    total = np.log(a / mu**2)
    for r in ((-b - sq) / (2 * a), (-b + sq) / (2 * a)):
        total += (_xlogx(1 - r) - _xlogx(-r) - 1).real      # r ln(-r) = -(-r) ln(-r)
    return float(-total)


def g3_at_mt(Mt: float, alpha_s: float) -> float:
    """MS-bar g3(mu = M_t) in the six-flavour SM: Buttazzo et al. eq. (3.5) --
    four-loop QCD running from M_Z, two-loop electroweak, three-loop matching."""
    return (1.1666 + 0.00314 * (alpha_s - REF["alpha_s"]) / 0.0007
            - 0.00046 * (Mt - REF["Mt"]))


def lambda_1loop(Mh, Mt, MW, MZ, V, mu):
    r = (3 * Mt**2 * (Mh**2 - 4 * Mt**2) * B0(Mh, Mt, Mt, mu) + 3 * Mh**2 * A0(Mt, mu)
         + 0.25 * (Mh**4 - 4 * Mh**2 * MZ**2 + 12 * MZ**4) * B0(Mh, MZ, MZ, mu)
         + Mh**2 * (7 * MW**2 - 4 * MZ**2) / (2 * (MZ**2 - MW**2)) * A0(MZ, mu)
         + 0.5 * (Mh**4 - 4 * Mh**2 * MW**2 + 12 * MW**4) * B0(Mh, MW, MW, mu)
         - 3 * Mh**2 * MW**2 / (2 * (Mh**2 - MW**2)) * A0(Mh, mu)
         + Mh**2 / 2 * (-11 + 3 * Mh**2 / (Mh**2 - MW**2) - 3 * MW**2 / (MZ**2 - MW**2))
         * A0(MW, mu)
         + 9 / 4 * Mh**4 * B0(Mh, Mh, Mh, mu)
         + 0.25 * (Mh**4 + Mh**2 * (MZ**2 + 2 * MW**2 - 6 * Mt**2)
                   - 8 * (MZ**4 + 2 * MW**4)))
    return K * r / V**4


def yt_1loop(Mh, Mt, MW, MZ, V, mu, g3sq):
    ew = (-(Mh**2 - 4 * Mt**2) * B0(Mt, Mh, Mt, mu)
          + (Mt**2 * (80 * MW**2 * MZ**2 - 64 * MW**4 - 7 * MZ**4) + 40 * MW**2 * MZ**4
             - 32 * MW**4 * MZ**2 - 17 * MZ**6) / (9 * Mt**2 * MZ**2) * B0(Mt, Mt, MZ, mu)
          + (Mt**2 * MW**2 + Mt**4 - 2 * MW**4) / Mt**2 * B0(Mt, 0.0, MW, mu)
          + (3 * Mh**2 / (Mh**2 - MW**2) + 2 * MW**2 / Mt**2 + 3 * MW**2 / (MW**2 - MZ**2)
             - 10) * A0(MW, mu)
          + (3 * MW**2 / (MW**2 - Mh**2) + 1) * A0(Mh, mu)
          + (36 * Mt**2 * MZ**2 - 56 * MW**2 * MZ**2 + 64 * MW**4 - 17 * MZ**4)
          / (9 * Mt**2 * MZ**2) * A0(Mt, mu)
          + (3 * MW**2 / (MZ**2 - MW**2) + (32 * MW**4 - 40 * MW**2 * MZ**2 + 17 * MZ**4)
             / (9 * Mt**2 * MZ**2) - 3) * A0(MZ, mu)
          + Mh**2 / 2 - 3 * Mt**2 - 9 * MW**2 + 7 * MZ**2 / 18 + 64 * MW**4 / (9 * MZ**2))
    qcd = g3sq * (-8 * A0(Mt, mu) / Mt**2 - 8 / 3)
    return K * Mt / (np.sqrt(2) * V) * (ew / V**2 + qcd)


def g2_1loop(Mh, Mt, MW, MZ, V, mu):
    r = ((Mh**4 / (6 * MW**2) - 2 * Mh**2 / 3 + 2 * MW**2) * B0(MW, Mh, MW, mu)
         + (-Mt**4 / MW**2 - Mt**2 + 2 * MW**2) * B0(MW, 0.0, Mt, mu)
         + (-48 * MW**4 / MZ**2 + MZ**4 / MW**2 - 68 * MW**2 + 16 * MZ**2) / 6
         * B0(MW, MW, MZ, mu)
         + (Mh**2 * (9 / (Mh**2 - MW**2) + 1 / MW**2) + MZ**2 / MW**2
            + MW**2 * (9 / (MW**2 - MZ**2) + 48 / MZ**2) - 27) / 6 * A0(MW, mu)
         + (2 - Mh**2 * (Mh**2 + 8 * MW**2) / (6 * MW**2 * (Mh**2 - MW**2))) * A0(Mh, mu)
         + (Mt**2 / MW**2 + 1) * A0(Mt, mu)
         + (24 * MW**2 / MZ**2 - MZ**2 / MW**2 + 9 * MW**2 / (MZ**2 - MW**2) - 17) / 6
         * A0(MZ, mu)
         + (-3 * Mh**2 + 18 * Mt**2 + 288 * MW**4 / MZ**2 - 374 * MW**2 - 3 * MZ**2) / 36)
    return K * 2 * MW / V**3 * r


def gY_1loop(Mh, Mt, MW, MZ, V, mu, last_term_mass4: bool = True):
    """One-loop correction to gY. The paper prints the term 36 MW^2/(MZ^2-MW^2)
    in its last bracket, which is dimensionless among mass-squared terms; with
    ``last_term_mass4`` it is read as 36 MW^4/(MZ^2-MW^2). test_physics.py shows
    that only this reading reproduces the paper's own one-loop value (table 3)."""
    t36 = 36 * MW**4 / (MZ**2 - MW**2) if last_term_mass4 else 36 * MW**2 / (MZ**2 - MW**2)
    r = ((88 / 9 - 124 * MW**2 / (9 * MZ**2) + (Mh**2 + 34 * MW**2) / (6 * (MZ**2 - MW**2)))
         * A0(MZ, mu)
         + (Mh**2 - 4 * MW**2) / (2 * (Mh**2 - MW**2)) * A0(Mh, mu)
         + (-7 / 9 - Mt**2 / (MZ**2 - MW**2) + 64 * MW**2 / (9 * MZ**2)) * A0(Mt, mu)
         + (Mh**4 + 2 * MW**2 * (MW**2 - 15 * MZ**2) + 3 * Mh**2 * (2 * MW**2 + 7 * MZ**2))
         / (6 * (Mh**2 - MW**2) * (MW**2 - MZ**2)) * A0(MW, mu)
         - (Mt**4 + MW**2 * Mt**2 - 2 * MW**4) / (MW**2 - MZ**2) * B0(MW, 0.0, Mt, mu)
         - (Mh**4 - 4 * MZ**2 * Mh**2 + 12 * MZ**4) / (6 * (MW**2 - MZ**2)) * B0(MZ, Mh, MZ, mu)
         + (Mh**4 - 4 * MW**2 * Mh**2 + 12 * MW**4) / (6 * (MW**2 - MZ**2)) * B0(MW, Mh, MW, mu)
         + (MZ**6 - 48 * MW**6 - 68 * MZ**2 * MW**4 + 16 * MZ**4 * MW**2)
         / (6 * MZ**2 * (MW**2 - MZ**2)) * B0(MW, MW, MZ, mu)
         + (-23 * MW**2 + 7 * Mt**2 + 17 * MZ**2 - 64 * Mt**2 * MW**2 / MZ**2
            - 9 * MW**2 * (Mt**2 - MW**2) / (MZ**2 - MW**2)) / 9 * B0(MZ, Mt, Mt, mu)
         + (MZ**6 - 48 * MW**6 - 68 * MZ**2 * MW**4 + 16 * MZ**4 * MW**2)
         / (6 * MZ**2 * (MZ**2 - MW**2)) * B0(MZ, MW, MW, mu)
         + (576 * MW**4 / MZ**2 - 242 * MW**2 - 3 * Mh**2 + 257 * MZ**2 + t36
            + Mt**2 * (82 - 256 * MW**2 / MZ**2)) / 36)
    return K * 2 * np.sqrt(MZ**2 - MW**2) / V**3 * r


def two_loop_terms(Mh: float, Mt: float, g3sq: float) -> dict:
    """Buttazzo et al.'s numerical two-loop thresholds at mu = M_t (section 2);
    linear fits around their inputs, small corrections in their own right."""
    dh, dt = Mh - 125.0, Mt - 173.0
    k2 = K * K
    return {
        "lam": k2 * (g3sq * (-23.88 + 0.12 * dh - 0.64 * dt) + (-9.45 - 0.12 * dh - 0.21 * dt)),
        "yt": k2 * ((6.48 - 0.01 * dh + 0.18 * dt) + g3sq * (-7.53 + 0.09 * dh - 0.23 * dt)
                    + g3sq**2 * (-145.08 - 0.84 * dt)),
        "g2": k2 * ((2.25 + 0.01 * dh + 0.01 * dt) + g3sq * (3.00 + 0.01 * dt)),
        "gY": k2 * ((-7.55 - 0.01 * dh - 0.11 * dt) + g3sq * (-14.66 - 0.14 * dt)),
    }


def alpha_s_5(mu: float, alpha_s_mz: float, MZ: float) -> float:
    """alpha_s in the five-flavour theory, four-loop running from M_Z."""
    nf = 5
    b0 = 11 - 2 / 3 * nf
    b1 = 102 - 38 / 3 * nf
    b2 = 2857 / 2 - 5033 / 18 * nf + 325 / 54 * nf**2
    b3 = (149753 / 6 + 3564 * ZETA3 - (1078361 / 162 + 6508 / 27 * ZETA3) * nf
          + (50065 / 162 + 6472 / 81 * ZETA3) * nf**2 + 1093 / 729 * nf**3)
    a = alpha_s_mz / (4 * np.pi)
    s0, s1 = 2 * np.log(MZ), 2 * np.log(mu)
    n = 400
    h = (s1 - s0) / n

    def f(x):
        return -(b0 * x**2 + b1 * x**3 + b2 * x**4 + b3 * x**5)

    for _ in range(n):
        k1 = f(a)
        k2 = f(a + 0.5 * h * k1)
        k3 = f(a + 0.5 * h * k2)
        k4 = f(a + h * k3)
        a += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return 4 * np.pi * a


def mb_running(mu: float, mb_mb: float, alpha_s_mz: float, MZ: float) -> float:
    """MS-bar m_b(mu) from m_b(m_b), three-loop QCD mass anomalous dimension
    (five flavours) -- only feeds the tiny bottom Yukawa coupling."""
    nf = 5
    b0 = 11 - 2 / 3 * nf
    b1 = 102 - 38 / 3 * nf
    b2 = 2857 / 2 - 5033 / 18 * nf + 325 / 54 * nf**2
    g0, g1 = 4.0, 202 / 3 - 20 / 9 * nf
    g2 = 1249 + (-2216 / 27 - 160 / 3 * ZETA3) * nf - 140 / 81 * nf**2
    # d ln m / d a = gamma_m(a) / beta(a), integrated from a(m_b) to a(mu)
    a1 = alpha_s_5(mb_mb, alpha_s_mz, MZ) / (4 * np.pi)
    a2 = alpha_s_5(mu, alpha_s_mz, MZ) / (4 * np.pi)
    xs = np.linspace(a1, a2, 2001)
    integrand = (g0 * xs + g1 * xs**2 + g2 * xs**3) / (b0 * xs**2 + b1 * xs**3 + b2 * xs**4)
    return mb_mb * np.exp(np.trapezoid(integrand, xs))


@dataclass
class Inputs:
    """Measured inputs (GeV); defaults are filled from config by the caller."""
    Mh: float
    Mt: float
    MW: float
    MZ: float
    alpha_s: float
    GF: float
    mb_mb: float = 4.183          # PDG, MS-bar m_b(m_b)
    Mtau: float = 1.77693
    extra: dict = field(default_factory=dict)

    @property
    def V(self) -> float:
        return (np.sqrt(2) * self.GF) ** -0.5


def match(inp: Inputs, order: str = "NNLO") -> dict:
    """MS-bar couplings at mu = M_t from the measured masses.

    order: "LO" (tree level), "NLO" (+ one loop), "NNLO" (+ two loops, and
    the three-loop pure-QCD term of y_t, as in Buttazzo et al.'s eq. 3.3).
    """
    Mh, Mt, MW, MZ, V = inp.Mh, inp.Mt, inp.MW, inp.MZ, inp.V
    mu = Mt
    g3 = g3_at_mt(Mt, inp.alpha_s)
    g3sq = g3 * g3
    lam = Mh**2 / (2 * V**2)
    yt = np.sqrt(2) * Mt / V
    g2 = 2 * MW / V
    gY = 2 * np.sqrt(MZ**2 - MW**2) / V
    if order in ("NLO", "NNLO"):
        lam += lambda_1loop(Mh, Mt, MW, MZ, V, mu)
        yt += yt_1loop(Mh, Mt, MW, MZ, V, mu, g3sq)
        g2 += g2_1loop(Mh, Mt, MW, MZ, V, mu)
        gY += gY_1loop(Mh, Mt, MW, MZ, V, mu)
    if order == "NNLO":
        t2 = two_loop_terms(Mh, Mt, g3sq)
        lam += t2["lam"]
        yt += t2["yt"]
        g2 += t2["g2"]
        gY += t2["gY"]
        a = g3sq / (4 * np.pi**2)                 # alpha_s^(6)(M_t) / pi
        yt += -80.4045 * a**3 * np.sqrt(2) * Mt / V   # three-loop pure QCD
    yb = np.sqrt(2) * mb_running(Mt, inp.mb_mb, inp.alpha_s, MZ) / V
    ytau = np.sqrt(2) * inp.Mtau / V
    return {"mu": mu, "lam": lam, "yt": yt, "g2": g2, "gY": gY, "g3": g3,
            "yb": yb, "ytau": ytau, "order": order}


def state(c: dict) -> np.ndarray:
    """Matching output -> RGE state vector (Gamma starts at 0 at mu = M_t)."""
    return np.array([5 / 3 * c["gY"] ** 2, c["g2"] ** 2, c["g3"] ** 2, c["yt"] ** 2,
                     c["yb"] ** 2, c["ytau"] ** 2, c["lam"], 0.0])


def butazzo_fit_formulas(Mh: float, Mt: float, alpha_s: float, MW: float) -> dict:
    """Buttazzo et al.'s published fits for the couplings at mu = M_t (sec. 3)."""
    da = (alpha_s - REF["alpha_s"]) / 0.0007
    dw = (MW - REF["MW"]) / REF["MW_err"]
    return {
        "lam": 0.12604 + 0.00206 * (Mh - REF["Mh"]) - 0.00004 * (Mt - REF["Mt"]),
        "yt": 0.93690 + 0.00556 * (Mt - REF["Mt"]) - 0.00042 * da,
        "g2": 0.64779 + 0.00004 * (Mt - REF["Mt"]) + 0.00011 * dw,
        "gY": 0.35830 + 0.00011 * (Mt - REF["Mt"]) - 0.00020 * dw,
        "g3": g3_at_mt(Mt, alpha_s),
    }


# ---------------------------------------------------------------------------
# The effective potential at large field values (appendix C)
# ---------------------------------------------------------------------------


def _li2(z):
    """Dilogarithm for complex arguments (scipy's spence(z) = Li2(1 - z))."""
    return spence(1 - np.asarray(z, dtype=complex))


def _xi(z):
    z = np.asarray(z, dtype=complex)
    s = np.sqrt(z * z - 4 * z)
    w = (z - s) / (2 * z)
    return s * (2 * np.log(w) ** 2 - np.log(z) ** 2 - 4 * _li2(w) + np.pi**2 / 3)


def lambda_eff(X: np.ndarray, loops: int = 2) -> np.ndarray:
    """lambda_eff(h) with V_eff = lambda_eff h^4 / 4 for h >> v, couplings at
    mu = h (Landau gauge; the two-loop part in the paper's lambda -> 0 limit,
    "well justified around the instability region")."""
    G1, G2, G3, YT, YB, YTA, L, Gam = X[:8]
    GY = 0.6 * G1
    rW = np.log(G2 / 4) + 2 * Gam
    rZ = np.log((G2 + GY) / 4) + 2 * Gam
    rt = np.log(YT / 2) + 2 * Gam
    one = K * (3 / 8 * G2**2 * (rW - 5 / 6) + 3 / 16 * (G2 + GY) ** 2 * (rZ - 5 / 6)
               - 3 * YT**2 * (rt - 3 / 2)
               + 3 * L**2 * (4 * np.log(np.abs(L)) - 6 + 3 * np.log(3) + 8 * Gam))
    total = L + (one if loops >= 1 else 0.0)
    if loops >= 2:
        rtW = (rt - rW) * (np.log((YT / 2 - G2 / 4).astype(complex)) + 2 * Gam)
        pi2 = np.pi**2
        two = (8 * G3 * YT**2 * (3 * rt**2 - 8 * rt + 9)
               + 0.5 * YT**3 * (-6 * rt * rW - 3 * rt**2 + 48 * rt - 6 * rtW - 69 - pi2)
               + 3 * YT * G2**2 / 16 * (8 * rW + 4 * rZ - 3 * rt**2 - 6 * rt * rZ - 12 * rt
                                        + 12 * rtW + 15 + 2 * pi2)
               + YT * GY**2 / 48 * (27 * rt**2 - 54 * rt * rZ - 68 * rt - 28 * rZ + 189)
               + YT * G2 * GY / 8 * (9 * rt**2 - 18 * rt * rZ + 4 * rt + 44 * rZ - 57)
               + G2**3 / 192 * (36 * rt * rZ + 54 * rt**2 - 414 * rW * rZ + 69 * rW**2
                                + 1264 * rW + 156 * rZ**2 + 632 * rZ - 144 * rtW - 2067
                                + 90 * pi2)
               + G2**2 * GY / 192 * (12 * rt * rZ - 6 * rt**2 - 6 * rW * (53 * rZ + 50)
                                     + 213 * rW**2 + 4 * rZ * (57 * rZ - 91) + 817 + 46 * pi2)
               + G2 * GY**2 / 576 * (132 * rt * rZ - 66 * rt**2 + 306 * rW * rZ - 153 * rW**2
                                     - 36 * rW + 924 * rZ**2 - 4080 * rZ + 4359 + 218 * pi2)
               + GY**3 / 576 * (6 * rZ * (34 * rt + 3 * rW - 470) - 102 * rt**2 - 9 * rW**2
                                + 708 * rZ**2 + 2883 + 206 * pi2)
               + YT**2 / 6 * (4 * GY * (3 * rt**2 - 8 * rt + 9) - 9 * G2 * (rt - rW + 1))
               + 0.75 * (G2**3 - 3 * G2**2 * YT + 4 * YT**3) * _li2(G2 / (2 * YT))
               + YT / 48 * _xi((G2 + GY) / (2 * YT))
               * (9 * G2**2 - 6 * G2 * GY + 17 * GY**2
                  + 2 * YT * (7 * GY - 73 * G2 + 64 * G2**2 / (GY + G2)))
               + G2 / 64 * _xi((G2 + GY) / G2)
               * (18 * G2 * GY + GY**2 - 51 * G2**2 - 48 * G2**3 / (GY + G2)))
        total = total + K * K * np.real(two)
    return np.exp(4 * Gam) * total


# ---------------------------------------------------------------------------
# Vacuum decay (Buttazzo et al. section 6.3)
# ---------------------------------------------------------------------------


def bounce_action(lam: np.ndarray) -> np.ndarray:
    """Tree-level action of the Fubini bounce, S = 8 pi^2 / (3 |lambda|)."""
    return 8 * np.pi**2 / (3 * np.abs(lam))
