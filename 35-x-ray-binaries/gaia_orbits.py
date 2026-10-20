"""Dormant black holes from Gaia DR3 orbits: Thiele-Innes constants -> orbit ->
astrometric mass function -> mass of the dark companion.

Gaia DR3 publishes a photocentre orbit as four Thiele-Innes constants (A, B,
F, G, in mas) plus e, P and the periastron time, with the full correlation
matrix (``corr_vec``). The geometric elements follow from (Binnendijk 1960;
Halbwachs et al. 2023, A&A 674, A9):

    a0^2 (1 + cos^2 i) = A^2 + B^2 + F^2 + G^2,   a0^2 cos i = AG - BF,
    omega + Omega = atan2(B - F, A + G),          omega - Omega = atan2(-B - F, A - G).

For a dark companion the photocentre is the star itself, so a0 is the
star's orbit about the centre of mass and

    (a0 / parallax)^3 / (P / yr)^2 = M2^3 / (M1 + M2)^2,

the astrometric mass function (AU, years, solar masses). With M1 from the
star's spectrum, M2 follows. We sample the DR3 covariance (parameter order as
in the DPAC ``nsstools`` package) and push every draw through these formulas.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from fetch_data import GAIA_NSS_PATH
from masses import astrometric_mass_function, draw, solve_m2_astrometric

BASE = ["ra", "dec", "parallax", "pmra", "pmdec",
        "a_thiele_innes", "b_thiele_innes", "f_thiele_innes", "g_thiele_innes",
        "c_thiele_innes", "h_thiele_innes"]
TAIL = {
    "Orbital": ["eccentricity", "period", "t_periastron"],
    "AstroSpectroSB1": ["center_of_mass_velocity", "eccentricity", "period", "t_periastron"],
}


def load_nss() -> pd.DataFrame:
    df = pd.read_csv(GAIA_NSS_PATH)
    df.columns = [c.lower() for c in df.columns]  # the archive's case varies
    return df


def covariance(row: pd.Series) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Parameter names, values and covariance of one DR3 solution.

    The parameters present are those of the solution type whose ``_error`` is
    finite; ``corr_vec`` is the upper triangle of their correlation matrix,
    column by column (c01, c02, c12, c03, ...)."""
    names = [n for n in BASE + TAIL[row["nss_solution_type"]]
             if pd.notna(row.get(f"{n}_error"))]
    corr = np.array([float(x) for x in str(row["corr_vec"]).strip("()[]").split(",")])
    n = len(names)
    if len(corr) != n * (n - 1) // 2:
        raise ValueError(f"{row['source_id']}: {len(corr)} correlations for {n} parameters")
    err = np.array([row[f"{k}_error"] for k in names], float)
    cmat = np.eye(n)
    k = 0
    for j in range(n):
        for i in range(j):
            cmat[i, j] = cmat[j, i] = corr[k]
            k += 1
    mean = np.array([row[k] for k in names], float)
    return names, mean, cmat * np.outer(err, err)


def campbell(A, B, F, G, H=None):
    """Thiele-Innes -> (a0 [mas], i [deg], omega [deg], Omega [deg]); vectorised.
    The astrometric orbit alone cannot tell (omega, Omega) from
    (omega + 180, Omega + 180): Omega is put in [0, 180) unless the
    radial-velocity constant H is given (AstroSpectroSB1 solutions), which
    fixes the sign as in DPAC's nsstools (flip when H cos(omega) < 0)."""
    A, B, F, G = (np.asarray(x, float) for x in (A, B, F, G))
    u = 0.5 * (A**2 + B**2 + F**2 + G**2)
    v = A * G - B * F
    a0 = np.sqrt(u + np.sqrt(np.clip(u**2 - v**2, 0, None)))
    incl = np.degrees(np.arccos(np.clip(v / a0**2, -1, 1)))
    s = np.arctan2(B - F, A + G)
    d = np.arctan2(-B - F, A - G)
    omega, Omega = 0.5 * (s + d), 0.5 * (s - d)
    flip = Omega < 0
    Omega = np.where(flip, Omega + np.pi, Omega)
    omega = np.where(flip, omega + np.pi, omega)
    if H is not None:
        flip = np.asarray(H) * np.cos(omega) < 0
        Omega = np.where(flip, Omega + np.pi, Omega)
        omega = np.where(flip, omega + np.pi, omega)
    omega = (omega + np.pi) % (2 * np.pi) - np.pi
    return a0, incl, np.degrees(omega), np.degrees(Omega)


def thiele_innes(a0, incl_deg, omega_deg, Omega_deg):
    """Campbell -> Thiele-Innes (A, B, F, G), the inverse of ``campbell``."""
    i, w, W = np.radians(incl_deg), np.radians(omega_deg), np.radians(Omega_deg)
    A = a0 * (np.cos(w) * np.cos(W) - np.sin(w) * np.sin(W) * np.cos(i))
    B = a0 * (np.cos(w) * np.sin(W) + np.sin(w) * np.cos(W) * np.cos(i))
    F = -a0 * (np.sin(w) * np.cos(W) + np.cos(w) * np.sin(W) * np.cos(i))
    G = -a0 * (np.sin(w) * np.sin(W) - np.cos(w) * np.cos(W) * np.cos(i))
    return A, B, F, G


def orbit_on_sky(A, B, F, G, e, n: int = 400):
    """Photocentre track over one orbit: (Delta alpha*, Delta delta) in mas.
    Delta alpha* = B X + G Y, Delta delta = A X + F Y, with X = cos E - e and
    Y = sqrt(1 - e^2) sin E (E the eccentric anomaly)."""
    E = np.linspace(0, 2 * np.pi, n)
    X, Y = np.cos(E) - e, np.sqrt(1 - e**2) * np.sin(E)
    return B * X + G * Y, A * X + F * Y


def analyse(name: str, source_id: int, m_star_spec, rng: np.random.Generator,
            n: int = 200_000) -> dict:
    """Monte Carlo over the DR3 covariance (and the star's mass) for one system."""
    df = load_nss()
    row = df[df["source_id"] == source_id].iloc[0]
    names, mean, cov = covariance(row)
    x = rng.multivariate_normal(mean, cov, size=n, method="eigh")
    col = {k: x[:, j] for j, k in enumerate(names)}
    a0, incl, omega, Omega = campbell(col["a_thiele_innes"], col["b_thiele_innes"],
                                      col["f_thiele_innes"], col["g_thiele_innes"],
                                      col.get("h_thiele_innes"))
    good = (col["parallax"] > 0) & (col["eccentricity"] >= 0) & (col["eccentricity"] < 1) \
        & (col["period"] > 0)
    f_ast = astrometric_mass_function(a0, col["parallax"], col["period"])
    m1 = draw(m_star_spec, n, rng, lo=0.0)
    m2 = solve_m2_astrometric(f_ast, m1)
    a0_c, incl_c, omega_c, Omega_c = campbell(
        row["a_thiele_innes"], row["b_thiele_innes"], row["f_thiele_innes"],
        row["g_thiele_innes"], row["h_thiele_innes"] if "h_thiele_innes" in names else None)

    def stats(v):
        v = v[good]
        lo, med, hi = np.quantile(v, [0.16, 0.5, 0.84])
        return {"median": float(med), "plus": float(hi - med), "minus": float(med - lo),
                "mean": float(v.mean()), "std": float(v.std())}

    return {
        "name": name, "source_id": int(source_id), "solution_type": row["nss_solution_type"],
        "parameters": names,
        "catalogue": {k: float(row[k]) for k in names + ["significance", "goodness_of_fit"]
                      if k in row and pd.notna(row[k])},
        "central": {"a0_mas": float(a0_c), "incl_deg": float(incl_c),
                    "omega_deg": float(omega_c), "Omega_deg": float(Omega_c)},
        "a0_mas": stats(a0), "incl_deg": stats(incl), "omega_deg": stats(omega),
        "Omega_deg": stats(Omega), "f_ast": stats(f_ast), "m2": stats(m2),
        "a_star_au": stats(a0 / col["parallax"]),
        "fraction_valid": float(good.mean()),
        "_draws": {"a0": a0[good][:4000], "incl": incl[good][:4000], "m2": m2[good],
                   "f_ast": f_ast[good],
                   "A": col["a_thiele_innes"][good][:300], "B": col["b_thiele_innes"][good][:300],
                   "F": col["f_thiele_innes"][good][:300], "G": col["g_thiele_innes"][good][:300],
                   "e": col["eccentricity"][good][:300]},
        "_row": row,
    }
