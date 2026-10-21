"""Sirius B as a white dwarf: its radius, and the mass-radius relation it tests.

* **Chandrasekhar's relation**, computed here from first principles: a star
  held up only by the pressure of a cold, fully degenerate electron gas (the
  exact special-relativistic equation of state), integrated outward from the
  centre for a range of central densities. With mu_e = 2 (helium, carbon or
  oxygen: two nucleons per electron) the mass tends to ~1.46 Msun as the
  radius shrinks to zero: the Chandrasekhar limit. Iron (mu_e = 2.15) gives
  smaller stars.
* **Realistic models**: the Montreal evolutionary sequences (Bedard et al.
  2020; C/O core, thick hydrogen envelope) interpolated to Sirius B's
  temperature. A 25,000 K star is ~1-2% larger than a cold one of the same mass.
* **Sirius B's radius from its brightness and temperature**: the Montreal
  model atmospheres give the surface brightness of a 25,000 K hydrogen
  atmosphere, so the measured V magnitude and the parallax fix the angular,
  then the linear, size. L = 4 pi R^2 sigma T^4 follows.
* **Gravitational redshift**: v = GM/(Rc) from our mass and radius, against
  Joyce et al.'s (2018) measured 80.65 +- 0.77 km/s (project N3's showcase).
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

import config as C
import literature as L

# cgs constants for the structure integration
ME, CL, HP, MU, GG = 9.1093837015e-28, 2.99792458e10, 6.62607015e-27, 1.66053906660e-24, 6.67430e-8
MSUN_G, RSUN_CM = C.M_SUN * 1e3, C.R_SUN * 1e2
K_P = np.pi * ME**4 * CL**5 / (3.0 * HP**3)          # pressure scale, dyn/cm^2
K_RHO = (8.0 * np.pi / 3.0) * (ME * CL / HP) ** 3 * MU   # density per unit mu_e at x = 1, g/cm^3


def _pressure(x):
    return K_P * (x * (2 * x * x - 3) * np.sqrt(1 + x * x) + 3 * np.arcsinh(x))


def chandrasekhar_star(xc, mu_e=2.0):
    """(M [Msun], R [Rsun]) of a cold degenerate star with central Fermi momentum x_c = p_F / (m_e c)."""
    rho_of = lambda x: mu_e * K_RHO * x**3
    dPdx = lambda x: 8.0 * K_P * x**4 / np.sqrt(1 + x * x)

    def rhs(r, y):
        x, m = y
        x = max(x, 0.0)
        rho = rho_of(x)
        return [-GG * m * rho / (r * r * dPdx(x)) if x > 0 else 0.0, 4 * np.pi * r * r * rho]

    def surface(r, y):
        return y[0] - 1e-5 * xc
    surface.terminal, surface.direction = True, -1
    # start slightly off the centre with the series solution
    rho_c = rho_of(xc)
    r0 = 1e-6 * (_pressure(xc) / (GG * rho_c**2)) ** 0.5
    m0 = 4.0 / 3.0 * np.pi * r0**3 * rho_c
    sol = solve_ivp(rhs, (r0, 1e11), [xc, m0], events=surface, rtol=1e-9, atol=[1e-12, 1e18],
                    max_step=np.inf, method="LSODA")
    R = sol.t_events[0][0] if len(sol.t_events[0]) else sol.t[-1]
    M = sol.y_events[0][0][1] if len(sol.t_events[0]) else sol.y[1, -1]
    return M / MSUN_G, R / RSUN_CM


def chandrasekhar_curve(mu_e=2.0, xcs=np.geomspace(0.05, 60.0, 70)):
    MR = np.array([chandrasekhar_star(x, mu_e) for x in xcs])
    return pd.DataFrame({"xc": xcs, "M": MR[:, 0], "R": MR[:, 1]})


# ---------------------------------------------------------------------------
# Montreal models
# ---------------------------------------------------------------------------


def read_sequence(path) -> pd.DataFrame:
    """One Bedard et al. 2020 cooling sequence: Teff, log g, R (cm), age (yr), L (erg/s)."""
    rows = []
    for ln in path.read_text().splitlines():
        parts = ln.split()
        if len(parts) == 6 and re.fullmatch(r"\d+", parts[0]):
            rows.append([float(p) for p in parts[1:]])
    return pd.DataFrame(rows, columns=["Teff", "logg", "R_cm", "age", "L"])


def montreal_mr(teff, masses=C.MONTREAL_MASSES):
    """Radius (Rsun) at a given Teff along each thick-H sequence -> R(M) at that temperature."""
    out = []
    for m in masses:
        s = read_sequence(C.MONTREAL_DIR / f"seq_{m:03d}_thick.txt")
        s = s.sort_values("Teff")
        if not (s["Teff"].min() <= teff <= s["Teff"].max()):
            continue
        R = np.interp(np.log10(teff), np.log10(s["Teff"]), s["R_cm"]) / RSUN_CM
        age = np.interp(np.log10(teff), np.log10(s["Teff"]), s["age"])
        out.append((m / 100.0, R, age))
    return pd.DataFrame(out, columns=["M", "R", "age"])


def table_da() -> pd.DataFrame:
    """Montreal synthetic photometry of pure-H atmospheres (Teff, log g, M, Mbol, V, Gaia G...)."""
    path = C.MONTREAL_DIR / "Table_DA"
    rows = []
    for ln in path.read_text().splitlines()[2:]:
        p = ln.split()
        if len(p) == 44:
            rows.append([float(v) for v in p])
    cols = ["Teff", "logg", "M", "Mbol", "BC", "U", "B", "V", "R", "I"]
    df = pd.DataFrame(np.array(rows)[:, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 38, 39, 40, 43]],
                      columns=cols + ["G3", "G3_BP", "G3_RP", "age"])
    g = 10 ** df["logg"] / 100.0                                   # m/s^2
    df["R_model"] = np.sqrt(C.G_NEWTON * df["M"] * C.M_SUN / g) / C.R_SUN
    return df


def surface_brightness(teff, logg, band="V"):
    """M_band + 5 log10(R/Rsun) of a DA atmosphere: depends on Teff (and weakly on log g)."""
    t = table_da()
    t["S"] = t[band] + 5 * np.log10(t["R_model"])
    out = []
    for lg in sorted(t["logg"].unique()):
        s = t[t["logg"] == lg].sort_values("Teff")
        out.append((lg, np.interp(np.log10(teff), np.log10(s["Teff"]), s["S"])))
    lgs, S = np.array(out).T
    return float(np.interp(logg, lgs, S))


def photometric_radius(mag, mag_err, plx, plx_err, teff, teff_err, logg, band="V", n=20000, seed=C.SEED):
    """Radius (Rsun) from apparent magnitude + parallax + model surface brightness (Monte Carlo)."""
    rng = np.random.default_rng(seed)
    m = rng.normal(mag, mag_err, n)
    p = rng.normal(plx, plx_err, n)
    t = rng.normal(teff, teff_err, n)
    # surface brightness is smooth in Teff: linearise around the nominal value
    S0 = surface_brightness(teff, logg, band)
    dS = (surface_brightness(teff + 200, logg, band) - surface_brightness(teff - 200, logg, band)) / 400.0
    S = S0 + dS * (t - teff)
    absmag = m + 5 * np.log10(p) + 5          # p in arcsec
    R = 10 ** (0.2 * (S - absmag))
    Lum = (R * C.R_SUN) ** 2 * 4 * np.pi * C.SIGMA_SB * t**4 / C.L_SUN
    return {"R": R, "L": Lum, "S0": S0, "dS_dT": dS}


def grav_redshift(M, R):
    """GM/(Rc) in km/s for M in Msun and R in Rsun."""
    return C.GM_SUN * np.asarray(M) / (np.asarray(R) * C.R_SUN * C.C_LIGHT) / 1e3


def parsons_table() -> pd.DataFrame:
    """Parsons et al. 2017 (MNRAS 470, 4473) Table: eclipsing white dwarfs with M and R."""
    tex = (C.LIT_DIR / "parsons2017_src" / "wd_mr.tex").read_text(encoding="latin-1")
    block = tex[tex.index("label{tab:wdparas}"):]
    block = block[:block.index("end{tabular}")]
    rows = []
    for ln in block.splitlines():
        if "\\pm" not in ln or ln.strip().startswith("Object"):
            continue
        cells = [c.strip() for c in ln.replace("\\\\", "").split("&")]
        def pm(s):
            v = re.findall(r"([-\d.]+)\s*\\pm\s*([\d.]+)", s)
            return (float(v[0][0]), float(v[0][1])) if v else (np.nan, np.nan)
        name = re.sub(r"\\,|\$|\{|\}", " ", cells[0]).replace("  ", " ").strip()
        teff, M, R = pm(cells[3]), pm(cells[4]), pm(cells[5])
        rows.append({"name": name, "Teff": teff[0], "M": M[0], "e_M": M[1], "R": R[0], "e_R": R[1]})
    return pd.DataFrame(rows)


def visual_binary_wds() -> pd.DataFrame:
    rows = [{"name": k, "Teff": v["Teff"][0], "M": v["M"][0], "e_M": v["M"][1], "R": v["R"][0], "e_R": v["R"][1]}
            for k, v in L.VISUAL_BINARY_WDS.items()]
    return pd.DataFrame(rows)
