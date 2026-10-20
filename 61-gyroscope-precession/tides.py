"""Tidal changes of the Earth's gravity field (IERS Conventions 2010, chapter 6).

The Moon and the Sun deform the solid Earth and move the oceans; both change the
Earth's external gravity field, which we express as small time-varying changes of
the normalised Stokes coefficients, added to the static field:

* **Solid Earth tides** (section 6.2): Step 1 with nominal anelastic Love
  numbers k_nm from the instantaneous lunar and solar positions (eq. 6.6, 6.7);
  Step 2, the frequency-dependent corrections for the diurnal (K1 above all),
  long-period and semidiurnal bands (Tables 6.5a-c); Step 3, the permanent tide
  removed because GOCO06s is a "zero tide" model (eq. 6.13-6.14).
* **Ocean tides** (section 6.3): FES2004's main waves as geopotential
  coefficients (eq. 6.15), completed with the secondary waves of Table 6.7 by
  linear admittance interpolation (eq. 6.16). This is what carries the 18.6-year
  nodal modulation of K1 and O1.
* **Pole tides** (sections 6.4, 6.5): the solid-Earth and ocean response to the
  centrifugal effect of polar motion (degree 2 only).

For the orbits of LAGEOS and LARES 2 the most important single term is K1: the
solid + ocean K1 tide makes each satellite's node oscillate with the period of
its own nodal regression (~2.9 years), with an amplitude of hundreds of
milliarcseconds, far above frame dragging's 31 mas/yr.
"""

from __future__ import annotations

import functools

import erfa
import numpy as np

import config as C
from gravity import legendre_bar

MJD0 = 2_400_000.5

# --- IERS Table 6.3: nominal Love numbers (anelastic Earth) -------------------
K2M = (0.30190 + 0.0j, 0.29830 - 0.00144j, 0.30102 - 0.00130j)
K3M = (0.093, 0.093, 0.093, 0.094)
K2M_PLUS = (-0.00089, -0.00080, -0.00057)
#: permanent part of the degree-2 zonal tide (eq. 6.14): A0 * H0 * k20
DC20_PERM = 4.4228e-8 * (-0.31460) * 0.30190

# --- Table 6.5a (diurnal, m=1): Doodson number, ip, op (1e-12) ----------------
TABLE_65A = [
    ("125.755", -0.1, 0.0), ("127.555", -0.1, 0.0), ("135.645", -0.1, 0.0), ("135.655", -0.7, 0.1),
    ("137.455", -0.1, 0.0), ("145.545", -1.3, 0.1), ("145.555", -6.8, 0.6), ("147.555", 0.1, 0.0),
    ("153.655", 0.1, 0.0), ("155.445", 0.1, 0.0), ("155.455", 0.4, 0.0), ("155.655", 1.3, -0.1),
    ("155.665", 0.3, 0.0), ("157.455", 0.3, 0.0), ("157.465", 0.1, 0.0), ("162.556", -1.9, 0.1),
    ("163.545", 0.5, 0.0), ("163.555", -43.4, 2.9), ("164.554", 0.6, 0.0), ("164.556", 1.6, -0.1),
    ("165.345", 0.1, 0.0), ("165.535", 0.1, 0.0), ("165.545", -8.8, 0.5), ("165.555", 470.9, -30.2),
    ("165.565", 68.1, -4.6), ("165.575", -1.6, 0.1), ("166.455", 0.1, 0.0), ("166.544", -0.1, 0.0),
    ("166.554", -20.6, -0.3), ("166.556", 0.3, 0.0), ("166.564", -0.3, 0.0), ("167.355", -0.2, 0.0),
    ("167.365", -0.1, 0.0), ("167.555", -5.0, 0.3), ("167.565", 0.2, 0.0), ("168.554", -0.2, 0.0),
    ("173.655", -0.5, 0.0), ("173.665", -0.1, 0.0), ("175.445", 0.1, 0.0), ("175.455", -2.1, 0.1),
    ("175.465", -0.4, 0.0), ("183.555", -0.2, 0.0), ("185.355", -0.1, 0.0), ("185.555", -0.6, 0.0),
    ("185.565", -0.4, 0.0), ("185.575", -0.1, 0.0), ("195.455", -0.1, 0.0), ("195.465", -0.1, 0.0),
]
# --- Table 6.5b (long period, m=0): Doodson number, ip, op (1e-12) ------------
TABLE_65B = [
    ("055.565", 16.6, -6.7), ("055.575", -0.1, 0.1), ("056.554", -1.2, 0.8), ("057.555", -5.5, 4.3),
    ("057.565", 0.1, -0.1), ("058.554", -0.3, 0.2), ("063.655", -0.3, 0.7), ("065.445", 0.1, -0.2),
    ("065.455", -1.2, 3.7), ("065.465", 0.1, -0.2), ("065.655", 0.1, -0.2), ("073.555", 0.0, 0.6),
    ("075.355", 0.0, 0.3), ("075.555", 0.6, 6.3), ("075.565", 0.2, 2.6), ("075.575", 0.0, 0.2),
    ("083.655", 0.1, 0.2), ("085.455", 0.4, 1.1), ("085.465", 0.2, 0.5), ("093.555", 0.1, 0.2),
    ("095.355", 0.1, 0.1),
]
# --- Table 6.5c (semidiurnal, m=2): Doodson number, amplitude (1e-12) ---------
TABLE_65C = [("245.655", -0.3), ("255.555", -1.2)]

# --- Table 6.7: astronomical amplitudes H_f (m) and pivot waves for FES2004 ---
# (Doodson number, H_f, pivot 1, pivot 2); main waves have no pivots.
TABLE_67 = [
    ("055.565", 0.02793, None, None), ("055.575", -0.00027, None, None),
    ("056.554", -0.00492, None, None), ("057.555", -0.03100, None, None),
    ("058.554", -0.00181, "057.555", "065.455"), ("063.655", -0.00673, "057.555", "065.455"),
    ("065.445", 0.00231, "057.555", "065.455"), ("065.455", -0.03518, None, None),
    ("065.465", 0.00229, "065.455", "075.555"), ("065.555", -0.00375, "065.455", "075.555"),
    ("065.655", 0.00188, "065.455", "075.555"), ("073.555", -0.00583, "065.455", "075.555"),
    ("075.355", -0.00288, "065.455", "075.555"), ("075.555", -0.06663, None, None),
    ("075.565", -0.02762, "075.555", "085.455"), ("075.575", -0.00258, "075.555", "085.455"),
    ("083.655", -0.00242, "075.555", "085.455"), ("083.665", -0.00100, "075.555", "085.455"),
    ("085.455", -0.01276, None, None), ("085.465", -0.00529, "085.455", "093.555"),
    ("093.555", -0.00204, None, None), ("095.355", -0.00169, "085.455", "093.555"),
    ("117.655", -0.00194, "135.455", "145.555"), ("125.755", -0.00664, "135.655", "145.555"),
    ("127.555", -0.00802, "135.655", "145.555"), ("135.645", -0.00947, "135.655", "145.555"),
    ("135.655", -0.05020, None, None), ("137.445", -0.00180, "135.655", "145.555"),
    ("137.455", -0.00954, "135.655", "145.555"), ("145.545", -0.04946, "135.655", "145.555"),
    ("145.555", -0.26221, None, None), ("145.755", 0.00170, "145.555", "165.555"),
    ("147.555", 0.00343, "145.555", "165.555"), ("153.655", 0.00194, "145.555", "165.555"),
    ("155.455", 0.00741, "145.555", "165.555"), ("155.555", -0.00399, "145.555", "165.555"),
    ("155.655", 0.02062, "145.555", "165.555"), ("155.665", 0.00414, "145.555", "165.555"),
    ("157.455", 0.00394, "145.555", "165.555"), ("162.556", -0.00714, "145.555", "165.555"),
    ("163.555", -0.12203, None, None), ("165.545", -0.00730, "145.555", "165.555"),
    ("165.555", 0.36878, None, None), ("165.565", 0.05001, "145.555", "165.555"),
    ("166.554", 0.00293, "145.555", "165.555"), ("167.555", 0.00525, "145.555", "165.555"),
    ("173.655", 0.00395, "145.555", "165.555"), ("175.455", 0.02062, "145.555", "165.555"),
    ("175.465", 0.00409, "145.555", "165.555"), ("183.555", 0.00342, "145.555", "165.555"),
    ("185.355", 0.00169, "145.555", "165.555"), ("185.555", 0.01129, "145.555", "165.555"),
    ("185.565", 0.00723, "145.555", "165.555"), ("195.455", 0.00216, "145.555", "165.555"),
    ("225.855", 0.00180, "235.755", "245.655"), ("227.655", 0.00467, "235.755", "245.655"),
    ("235.755", 0.01601, None, None), ("237.555", 0.01932, "235.755", "245.655"),
    ("245.555", -0.00389, "237.755", "245.655"), ("245.645", -0.00451, "237.755", "245.655"),
    ("245.655", 0.12099, None, None), ("247.455", 0.02298, "245.655", "255.555"),
    ("253.755", -0.00190, "245.655", "255.555"), ("254.556", -0.00218, "245.655", "255.555"),
    ("255.545", -0.02358, "245.655", "255.555"), ("255.555", 0.63192, None, None),
    ("256.554", 0.00192, "255.555", "275.555"), ("263.655", -0.00466, "255.555", "275.555"),
    ("265.455", -0.01786, "255.555", "275.555"), ("265.555", 0.00359, "255.555", "275.555"),
    ("265.655", 0.00447, "255.555", "275.555"), ("265.665", 0.00197, "255.555", "275.555"),
    ("272.556", 0.01720, "255.555", "275.555"), ("273.555", 0.29400, None, None),
    ("274.554", -0.00246, "255.555", "275.555"), ("275.555", 0.07996, None, None),
    ("275.565", 0.02383, "255.555", "275.555"), ("275.575", 0.00259, "255.555", "275.555"),
    ("285.455", 0.00447, "255.555", "275.555"), ("285.465", 0.00195, "255.555", "275.555"),
    ("455.555", None, None, None),
]


#: Two pivot waves in Table 6.7 are not FES2004 main waves (typos in the table):
#: 135.455 for 117.655 must be Q1 (135.655), and 237.755 for the two lines next to
#: N2 must be 2N2 (235.755), the main wave below them.
PIVOT_FIX = {"135.455": "135.655", "237.755": "235.755"}


def doodson_multipliers(code: str) -> np.ndarray:
    """'165.555' -> (1, 1, 0, 0, 0, 0): multipliers of (tau, s, h, p, N', ps)."""
    digits = code.replace(".", "").replace(",", "")
    n = [int(digits[0])] + [int(d) - 5 for d in digits[1:]]
    return np.array(n, float)


def doodson_arguments(mjd_tt, mjd_ut1) -> np.ndarray:
    """Doodson's fundamental arguments (N,6) = (tau, s, h, p, N', ps) in radians.

    From the Delaunay arguments of nutation theory (IERS 2010 ch.5) and GMST:
    s = F + Om, h = s - D, p = s - l, N' = -Om, ps = s - D - l', tau = GMST + pi - s.
    """
    mjd_tt = np.atleast_1d(np.asarray(mjd_tt, float))
    mjd_ut1 = np.atleast_1d(np.asarray(mjd_ut1, float))
    t = (mjd_tt - 51_544.5) / 36_525.0
    l, lp, F, D, Om = erfa.fal03(t), erfa.falp03(t), erfa.faf03(t), erfa.fad03(t), erfa.faom03(t)
    gmst = erfa.gmst06(MJD0, mjd_ut1, MJD0, mjd_tt)
    s = F + Om
    h = s - D
    p = s - l
    Np = -Om
    ps = s - D - lp
    tau = gmst + np.pi - s
    return np.stack([tau, s, h, p, Np, ps], axis=1)


def doodson_rates() -> np.ndarray:
    """Rates (rad/day) of the Doodson arguments, by differencing over one day at J2025."""
    a0 = doodson_arguments(60_676.0, 60_676.0)[0]
    a1 = doodson_arguments(60_677.0, 60_677.0)[0]
    d = (a1 - a0 + np.pi) % (2 * np.pi) - np.pi
    d[0] += 2 * np.pi        # tau advances by ~ one turn per day
    return d


# ---------------------------------------------------------------------------
# Ocean tides: FES2004 + admittance waves
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=4)
def ocean_waves(nmax: int = C.OCEAN_TIDE_DEGREE):
    """All ocean-tide waves: multipliers (W,6) and coefficient arrays (W, nmax+1, nmax+1, 4).

    Coefficients are C+, S+, C-, S- in eq. 6.15 (normalised), main waves from the
    FES2004 file, secondary waves interpolated by eq. 6.16.
    """
    main = {}
    with open(C.TIDE_DIR / "fes2004_Cnm-Snm.dat") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 8:
                continue
            try:
                n, m = int(parts[2]), int(parts[3])
                vals = [float(v) * 1e-11 for v in parts[4:]]
            except ValueError:
                continue
            code = f"{float(parts[0]):07.3f}"      # ' 55.565' -> '055.565'
            if n > nmax:
                continue
            arr = main.setdefault(code, np.zeros((nmax + 1, nmax + 1, 4)))
            arr[n, m] = vals
    rates = doodson_rates()
    H = {code: h for code, h, _, _ in TABLE_67}
    waves = list(main.items())
    for code, h, p1, p2 in TABLE_67:
        if p1 is None or code in main:
            continue
        p1, p2 = PIVOT_FIX.get(p1, p1), PIVOT_FIX.get(p2, p2)
        f = doodson_multipliers(code) @ rates
        f1 = doodson_multipliers(p1) @ rates
        f2 = doodson_multipliers(p2) @ rates
        arr = ((f - f1) / (f2 - f1) * (h / H[p2]) * main[p2]
               + (f2 - f) / (f2 - f1) * (h / H[p1]) * main[p1])
        waves.append((code, arr))
    codes = [w[0] for w in waves]
    mult = np.array([doodson_multipliers(c) for c in codes])
    coef = np.array([w[1] for w in waves])
    return codes, mult, coef


#: The 18.6-year (and 9.3-year) equilibrium ocean tides as given in the IERS file
#: come out in anti-phase with the solid-Earth response to the same potential
#: (C20 = -6.6e-11 cos N' against the solid tide's +3.9e-10 cos N'). An equilibrium
#: ocean must rise and fall with the solid tide; the IERS text itself notes an
#: "ongoing discussion" about this term's pi/2 phase (2011 update). With
#: om1_flip=True the two waves are put in phase with the solid tide.
NODAL_WAVES = ("055.565", "055.575")


def ocean_tide_coefficients(beta: np.ndarray, nmax: int = C.OCEAN_TIDE_DEGREE, om1_flip: bool = False):
    """Delta Cbar, Delta Sbar (T, nmax+1, nmax+1) from eq. 6.15 at Doodson arguments beta (T,6)."""
    codes, mult, coef = ocean_waves(nmax)
    if om1_flip:
        coef = coef.copy()
        for k, code in enumerate(codes):
            if code in NODAL_WAVES:
                coef[k] = -coef[k]
    theta = beta @ mult.T                                 # (T, W)
    c, s = np.cos(theta), np.sin(theta)
    Cp, Sp, Cm, Sm = (coef[..., k] for k in range(4))    # (W, n, m)
    dC = np.einsum("tw,wnm->tnm", c, Cp + Cm) + np.einsum("tw,wnm->tnm", s, Sp + Sm)
    dS = np.einsum("tw,wnm->tnm", c, Sp - Sm) - np.einsum("tw,wnm->tnm", s, Cp - Cm)
    dS[:, :, 0] = 0.0                                     # FES convention for zonal terms
    return dC, dS


# ---------------------------------------------------------------------------
# Solid Earth tides
# ---------------------------------------------------------------------------

def solid_tide_coefficients(sun_itrs: np.ndarray, moon_itrs: np.ndarray, beta: np.ndarray,
                            zero_tide: bool = True):
    """Delta Cbar, Delta Sbar (T, 5, 5) for n <= 4 from IERS Steps 1-3.

    sun_itrs, moon_itrs: (T,3) Earth-fixed positions (m); beta: Doodson arguments (T,6).
    """
    T = sun_itrs.shape[0]
    dC = np.zeros((T, 5, 5))
    dS = np.zeros((T, 5, 5))
    for pos, gm in ((moon_itrs, C.GM_MOON), (sun_itrs, C.GM_SUN)):
        r = np.linalg.norm(pos, axis=1)
        sinphi = pos[:, 2] / r
        lam = np.arctan2(pos[:, 1], pos[:, 0])
        for t in range(T):
            P = legendre_bar(3, sinphi[t])
            ratio = gm / C.GM_EARTH
            for n in (2, 3):
                fac = ratio * (C.R_EARTH / r[t]) ** (n + 1) / (2 * n + 1)
                for m in range(n + 1):
                    k = K2M[m] if n == 2 else K3M[m]
                    k = complex(k)
                    base = fac * P[n, m]
                    cm, sm = np.cos(m * lam[t]), np.sin(m * lam[t])
                    dC[t, n, m] += base * (k.real * cm + k.imag * sm)
                    dS[t, n, m] += base * (k.real * sm - k.imag * cm)
            fac4 = ratio * (C.R_EARTH / r[t]) ** 3 / 5.0
            for m in range(3):
                cm, sm = np.cos(m * lam[t]), np.sin(m * lam[t])
                dC[t, 4, m] += fac4 * K2M_PLUS[m] * P[2, m] * cm
                dS[t, 4, m] += fac4 * K2M_PLUS[m] * P[2, m] * sm
    # Step 2: frequency dependence
    for code, ip, op in TABLE_65B:
        th = beta @ doodson_multipliers(code)
        dC[:, 2, 0] += (ip * np.cos(th) - op * np.sin(th)) * 1e-12
    for code, ip, op in TABLE_65A:
        th = beta @ doodson_multipliers(code)
        dC[:, 2, 1] += (ip * np.sin(th) + op * np.cos(th)) * 1e-12
        dS[:, 2, 1] += (ip * np.cos(th) - op * np.sin(th)) * 1e-12
    for code, amp in TABLE_65C:
        th = beta @ doodson_multipliers(code)
        dC[:, 2, 2] += amp * np.cos(th) * 1e-12
        dS[:, 2, 2] += -amp * np.sin(th) * 1e-12
    # Step 3: zero-tide field already contains the permanent tide
    if zero_tide:
        dC[:, 2, 0] -= DC20_PERM
    dS[:, :, 0] = 0.0
    return dC, dS


def mean_pole(mjd) -> tuple[np.ndarray, np.ndarray]:
    """IERS secular mean pole (2018 update of the Conventions), arcsec."""
    t = (np.asarray(mjd, float) - 51_544.5) / 365.25
    return 0.0550 + 0.001677 * t, 0.3205 + 0.003460 * t


def pole_tide_coefficients(mjd, xp_rad, yp_rad):
    """Delta Cbar21, Delta Sbar21 from the solid-Earth (6.4) and ocean (6.5) pole tides."""
    xbar, ybar = mean_pole(mjd)
    as_ = 180 / np.pi * 3600
    m1 = xp_rad * as_ - xbar
    m2 = -(yp_rad * as_ - ybar)
    dC = -1.333e-9 * (m1 + 0.0115 * m2) - 2.1778e-10 * (m1 - 0.01724 * m2)
    dS = -1.333e-9 * (m2 - 0.0115 * m1) - 1.7232e-10 * (m2 - 0.03365 * m1)
    return dC, dS
