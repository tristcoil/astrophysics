"""Stellar masses from absolute G magnitudes.

The wide binary test needs each star's mass, and in a test of gravity the
masses cannot come from the orbits themselves, so they come from luminosity:
a main-sequence star's brightness fixes its mass. The relations used here:

* ``banik``: Banik et al. (2024, eq. 4): a cubic fit to M_G against
  x = ln(M/Msun), built from the Pecaut & Mamajek (2013) dwarf sequence (M_V
  and V - I_C) and the Riello et al. (2021) G - V transformation; used for
  0.6 < M_G < 11.1 with straight-line extensions of slope -3.3 (faint end)
  and -2.9 (bright end). Optionally with their small FLAME calibration (eq. 7).
* ``chae``: Chae (2023) table 1, log10 M = sum a_i M_G^i (the V-band based
  relation, his standard choice).
* ``ps23``: Pittordis & Sutherland (2023) eq. 3, M = 10^(0.074 (4.69 - M_G)).

We also rebuild the underlying (M_G, M) points from the Pecaut & Mamajek
table so the fits can be checked against what they summarise.
"""

from __future__ import annotations

import numpy as np

from config import RAW_DIR

PM_TABLE = RAW_DIR / "pecaut_mamajek" / "EEM_dwarf_UBVIJHK_colors_Teff.txt"

# --------------------------------------------------------------------------
# Banik et al. (2024)
# --------------------------------------------------------------------------
_C = (4.887, -5.693, 0.4164, 0.9611)
MG_FAINT, MG_BRIGHT = 11.1, 0.6
SLOPE_FAINT, SLOPE_BRIGHT = -3.3, -2.9


def _cubic(x):
    return _C[0] + _C[1] * x + _C[2] * x**2 + _C[3] * x**3


#: The cubic's turning points: it peaks at M_G = 11.13 (x = -1.557), which is
#: why Banik et al. switch to a straight line at M_G = 11.1; the other turning
#: point (x = +1.268) lies beyond the bright-end switch at M_G = 0.6.
X_TURN_FAINT = (-2 * _C[2] - (4 * _C[2]**2 - 12 * _C[3] * _C[1]) ** 0.5) / (6 * _C[3])
X_TURN_BRIGHT = (-2 * _C[2] + (4 * _C[2]**2 - 12 * _C[3] * _C[1]) ** 0.5) / (6 * _C[3])


def _cubic_root(MG: float) -> float:
    """x with cubic(x) = MG on the monotonic branch between the turning points."""
    lo, hi = X_TURN_FAINT, X_TURN_BRIGHT
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _cubic(mid) > MG:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


X_FAINT = _cubic_root(MG_FAINT)      # ln M where the cubic reaches M_G = 11.1
X_BRIGHT = _cubic_root(MG_BRIGHT)    # ... and M_G = 0.6


def banik_MG(mass):
    """Absolute G magnitude of a main-sequence star of ``mass`` (Msun)."""
    x = np.log(np.asarray(mass, float))
    out = _cubic(x)
    out = np.where(x < X_FAINT, MG_FAINT + SLOPE_FAINT * (x - X_FAINT), out)
    out = np.where(x > X_BRIGHT, MG_BRIGHT + SLOPE_BRIGHT * (x - X_BRIGHT), out)
    return out


def banik_mass(MG, flame: bool = True):
    """Inverse of banik_MG (bisection on the monotonic relation), optionally
    with Banik et al.'s FLAME calibration of eq. (7)."""
    MG = np.asarray(MG, float)
    lo = np.full(MG.shape, -6.0)
    hi = np.full(MG.shape, 4.0)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        brighter = banik_MG(np.exp(mid)) < MG     # mid too massive
        hi = np.where(brighter, mid, hi)
        lo = np.where(brighter, lo, mid)
    m = np.exp(0.5 * (lo + hi))
    if flame:
        m = flame_adjust(m)
    return m


def flame_adjust(m_pm):
    """Banik et al. (2024) eq. (7): M_F = M_PM - 0.07 tanh((M_PM - 0.75)/0.16) f."""
    m = np.asarray(m_pm, float)
    f = np.where(m <= 0.5, np.exp(2.0 - 1.0 / m),
                 np.where(m <= 0.75, 1.0, np.exp((0.75 - m) / 2.0)))
    return m - 0.07 * np.tanh((m - 0.75) / 0.16) * f


# --------------------------------------------------------------------------
# Chae (2023), table 1 (V-band based M_G); valid for his 4 < M_G < 14 range
# --------------------------------------------------------------------------
CHAE_V = (5.2951695081428651e-01, -1.5827136745981818e-01, 8.4871478522417707e-03,
          7.8449380571379954e-04, -5.2267549153639953e-05, -1.6957228195696253e-05,
          1.6858627515537989e-06, -6.8083605428022648e-08, 1.4781005839376326e-09,
          7.7057036188153745e-11, -4.4776519490406922e-12)
CHAE_J = (4.5004396006609515e-01, -7.5227632902604175e-02, -1.6733691959840702e-02,
          3.2486543639338823e-03, 9.7481038895336188e-06, -3.7585795718404064e-05,
          2.1882376017987368e-06, -2.6298611913795649e-08, 1.3151945936755533e-09,
          -9.3661053655257917e-11, 4.9792442749063264e-13)


def chae_mass(MG, which: str = "V"):
    coef = CHAE_V if which == "V" else CHAE_J
    MG = np.asarray(MG, float)
    return 10.0 ** np.polynomial.polynomial.polyval(MG, coef)


def ps23_mass(MG):
    """Pittordis & Sutherland (2023) eq. (3)."""
    return 10.0 ** (0.074 * (4.69 - np.asarray(MG, float)))


def mass_from_MG(MG, relation: str = "banik"):
    if relation == "banik":
        return banik_mass(MG, flame=True)
    if relation == "banik_noflame":
        return banik_mass(MG, flame=False)
    if relation == "chae":
        return chae_mass(MG, "V")
    if relation == "chae_J":
        return chae_mass(MG, "J")
    if relation == "ps23":
        return ps23_mass(MG)
    raise ValueError(relation)


# --------------------------------------------------------------------------
# Light from mass, for blending an unresolved companion's light into a star
# --------------------------------------------------------------------------

def g_luminosity(mass):
    """G-band luminosity (arbitrary units) of a main-sequence star."""
    return 10.0 ** (-0.4 * banik_MG(mass))


def blended_mass_table(n_m: int = 400, n_q: int = 201):
    """For a star seen at estimated mass M_est that is really an unresolved
    pair (primary M_p, companion q M_p), tabulate on a grid of (log M_est, q):
    M_p / M_est, the companion's light fraction L~, and the true/estimated
    mass ratio. 'Estimated' means what an observer infers from the blended
    light with the same relation (Banik et al. 2024, Section 3.2.1)."""
    logm = np.linspace(np.log10(0.06), np.log10(4.5), n_m)
    q = np.linspace(0.0, 1.0, n_q)
    Lm, Q = np.meshgrid(10**logm, q, indexing="ij")
    # bisection in M_p so that L(M_p) + L(q M_p) = L(M_est)
    target = g_luminosity(Lm)
    lo = np.log(Lm) - 1.0
    hi = np.log(Lm) + 1e-12
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        mp = np.exp(mid)
        lt = g_luminosity(mp) + np.where(Q > 0, g_luminosity(np.maximum(Q * mp, 1e-4)), 0.0)
        too_bright = lt > target
        hi = np.where(too_bright, mid, hi)
        lo = np.where(too_bright, lo, mid)
    mp = np.exp(0.5 * (lo + hi))
    Lc = np.where(Q > 0, g_luminosity(np.maximum(Q * mp, 1e-4)), 0.0)
    Lfrac = Lc / (g_luminosity(mp) + Lc)
    return dict(logm=logm, q=q, mp_over_mest=mp / Lm, L_frac=Lfrac,
                mtrue_over_mest=mp * (1 + Q) / Lm)


# --------------------------------------------------------------------------
# The Pecaut & Mamajek dwarf sequence itself
# --------------------------------------------------------------------------

def riello_G_minus_V(VI):
    """Riello et al. (2021) table C2, first Johnson-Cousins relation."""
    VI = np.asarray(VI, float)
    return -0.01597 - 0.02809 * VI - 0.2483 * VI**2 + 0.03656 * VI**3 - 0.002939 * VI**4


def pecaut_points():
    """(mass, M_G) of each dwarf spectral type with tabulated M_V, V - I_C and
    mass, M_G from M_V and Riello's G - V(V - I_C)."""
    rows = []
    started = False
    for line in PM_TABLE.read_text().splitlines():
        if line.startswith("#SpT"):
            if started:
                break
            started = True
            header = line[1:].split()
            continue
        if not started or line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < len(header) - 1:
            continue
        rec = dict(zip(header, parts))
        try:
            mv, vi, m = float(rec["Mv"]), float(rec["V-Ic"]), float(rec["Msun"])
        except (KeyError, ValueError):
            continue
        rows.append((rec["SpT"], m, mv, vi))
    spt = np.array([r[0] for r in rows])
    m = np.array([r[1] for r in rows])
    mv = np.array([r[2] for r in rows])
    vi = np.array([r[3] for r in rows])
    return dict(spt=spt, mass=m, MV=mv, VI=vi, MG=mv + riello_G_minus_V(vi))


if __name__ == "__main__":
    print("x_faint, x_bright:", X_FAINT, X_BRIGHT, np.exp(X_FAINT), np.exp(X_BRIGHT))
    for MG in (2.0, 4.0, 4.7, 6.0, 8.0, 10.0, 12.0, 14.0):
        print(MG, [round(float(mass_from_MG(MG, r)), 3) for r in ("banik", "banik_noflame", "chae", "ps23")])
    p = pecaut_points()
    sel = (p["MG"] > 0.6) & (p["MG"] < 11.1)
    print(len(p["mass"]), "PM points; rms of banik cubic in M_G over its range:",
          np.sqrt(np.mean((banik_MG(p["mass"][sel]) - p["MG"][sel]) ** 2)))
    t = blended_mass_table()
    i = np.argmin(abs(10**t["logm"] - 1.19))
    print("two 1 Msun stars appear as", 10**t["logm"][i], "-> true/est", t["mtrue_over_mest"][i, -1],
          "(Banik: 1.19 Msun, dM = 0.68)")
