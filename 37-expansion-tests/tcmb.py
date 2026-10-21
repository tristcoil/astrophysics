"""Test 2: was the cosmic microwave background hotter in the past?

If space expands, the CMB's photons are stretched along with it and the
radiation cools as T = T0 (1 + z). If redshifts came from light tiring on its
way to us, the universe would be static and the background would have the same
temperature everywhere, at all distances: T = T0.

Two kinds of thermometers read the CMB temperature at a distance:

* **Galaxy clusters (the Sunyaev-Zel'dovich effect).** Hot electrons in a
  cluster scatter some CMB photons to higher energies. Seen through the
  cluster, the CMB is dimmer below 217 GHz and brighter above. The crossover
  sits at x = h nu / (k T) = 3.83, where T is the CMB temperature *at the
  cluster*. That light then reaches us redshifted by (1 + z). If T = T0 (1 + z),
  the two factors cancel and every cluster shows its crossover at 217 GHz. If
  T = T0, a cluster at z = 0.5 would show it at 145 GHz.
* **Molecules and atoms in gas clouds in front of quasars.** CMB photons
  excite rotational levels of CO, HCN, ... and fine-structure levels of
  carbon. The ratio of populations is a thermometer: an excitation
  temperature. Collisions in the gas also excite the levels, so a measured
  excitation temperature is an upper limit unless corrected for them.

Here the published temperatures (``literature/tcmb_measurements.csv``) are
fitted with T(z) = T0 (1 + z)^(1 - beta): beta = 0 is expansion, beta = 1 is
tired light.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import LITERATURE_DIR, PAPER_DIR, T0_FIRAS

H_PLANCK = 6.62607015e-34
K_BOLTZ = 1.380649e-23
C_LIGHT = 2.99792458e8
#: Zero of the non-relativistic thermal SZ spectrum: x = h nu / kT = 3.830.
X_SZ_NULL = 3.8300


def load_measurements() -> pd.DataFrame:
    df = pd.read_csv(LITERATURE_DIR / "tcmb_measurements.csv", comment="#")
    df["sets"] = df["sets"].str.split(";")
    return df


def select(df, name):
    return df[df["sets"].apply(lambda s: name in s)].reset_index(drop=True)


def t_model(z, beta, t0=T0_FIRAS[0]):
    return t0 * (1.0 + np.asarray(z)) ** (1.0 - beta)


def _arrays(df):
    d = df[df["kind"] != "upper"]
    return (d["z"].to_numpy(float), d["T"].to_numpy(float), d["err_plus"].to_numpy(float),
            d["err_minus"].to_numpy(float), d["err_sys"].to_numpy(float))


def chi2_grid(df, betas, t0s):
    """Split-normal chi-square on a grid of (beta, T0) pairs (broadcast together).

    The upper error applies when the model is above the measurement, the lower
    error when below; quoted systematics add in quadrature. Upper limits are
    not used in fits (they are shown in figures).
    """
    z, T, ep, em, es = _arrays(df)
    betas = np.atleast_1d(np.asarray(betas, float))[:, None]
    t0s = np.atleast_1d(np.asarray(t0s, float))[:, None]
    model = t0s * (1.0 + z[None, :]) ** (1.0 - betas)
    resid = model - T[None, :]
    err = np.hypot(np.where(resid > 0, ep, em), es)
    return np.sum((resid / err) ** 2, axis=1)


def chi2(df, beta, t0=T0_FIRAS[0]):
    return float(chi2_grid(df, [beta], [t0])[0])


def _interval(grid, c):
    i = int(np.argmin(c))
    inside = grid[c <= c[i] + 1.0]
    return float(grid[i]), float(inside.max() - grid[i]), float(grid[i] - inside.min()), float(c[i])


def fit_beta(df, *, t0=T0_FIRAS[0], grid=np.linspace(-0.3, 0.3, 60001)):
    """Best beta and the Delta chi^2 = 1 interval (T0 fixed at FIRAS)."""
    best, plus, minus, c = _interval(grid, chi2_grid(df, grid, t0))
    return {"beta": best, "plus": plus, "minus": minus, "chi2": c,
            "dof": int((df["kind"] != "upper").sum() - 1)}


def fit_t0(df, beta=0.0):
    """T0 from the high-redshift data alone (FIRAS rows removed), beta fixed."""
    d = df[df["z"] > 0]
    grid = np.linspace(2.5, 3.0, 50001)
    best, plus, minus, c = _interval(grid, chi2_grid(d, beta, grid))
    return {"T0": best, "plus": plus, "minus": minus, "chi2": c,
            "n": int((d["kind"] != "upper").sum())}


def compare_models(df):
    """Chi-square of expansion (beta = 0) and tired light (beta = 1)."""
    n = int((df["kind"] != "upper").sum())
    return {"n": n, "chi2_expansion": chi2(df, 0.0), "chi2_tired": chi2(df, 1.0)}


# ---------------------------------------------------------------------------
# What the SZ measurements mean: where the cluster's shadow turns into a glow
# ---------------------------------------------------------------------------

def sz_spectrum(x):
    """Non-relativistic thermal SZ intensity change, in units of y * I0.

    Delta I / (y I0) = x^4 e^x / (e^x - 1)^2 * (x coth(x/2) - 4), with
    x = h nu / kT the frequency in units of the local CMB temperature.
    """
    x = np.asarray(x, float)
    ex = np.exp(x)
    return x ** 4 * ex / (ex - 1) ** 2 * (x / np.tanh(x / 2) - 4)


def null_frequency_ghz(t_local, z):
    """Observed frequency of the SZ null for a cluster whose CMB is at t_local."""
    return X_SZ_NULL * K_BOLTZ * np.asarray(t_local) / H_PLANCK / (1 + np.asarray(z)) / 1e9


def observed_null(df):
    """Each SZ measurement re-expressed as the observed crossover frequency."""
    d = df[df["method"] == "tSZ"].copy()
    d["nu_null"] = null_frequency_ghz(d["T"], d["z"])
    scale = d["nu_null"] / d["T"]
    d["nu_err_plus"] = np.hypot(d["err_plus"], d["err_sys"]) * scale
    d["nu_err_minus"] = np.hypot(d["err_minus"], d["err_sys"]) * scale
    return d


# ---------------------------------------------------------------------------
# Why molecular thermometers read warm: CO in statistical equilibrium
# ---------------------------------------------------------------------------

def read_lamda(path=PAPER_DIR / "lamda_co.dat"):
    """Parse a LAMDA molecular data file: levels, radiative and collisional rates."""
    lines = [l.strip() for l in path.read_text().splitlines()]
    it = iter(lines)

    def next_data():
        for l in it:
            if l and not l.startswith("!"):
                return l
        raise EOFError

    next_data(); next_data()  # molecule name, weight
    n_lev = int(next_data())
    lev = np.array([[float(v) for v in next_data().split()[:3]] for _ in range(n_lev)])
    energy_cm, g = lev[:, 1], lev[:, 2]
    n_rad = int(next_data())
    rad = np.array([[float(v) for v in next_data().split()[:5]] for _ in range(n_rad)])
    n_part = int(next_data())
    partners = {}
    for _ in range(n_part):
        name = next_data()
        n_coll = int(next_data())
        n_temp = int(next_data())
        temps = np.array([float(v) for v in next_data().split()])
        rows = np.array([[float(v) for v in next_data().split()] for _ in range(n_coll)])
        key = "pH2" if "pH2" in name or "para" in name.lower() else (
            "oH2" if "oH2" in name or "ortho" in name.lower() else name)
        partners[key] = {"temps": temps, "up": rows[:, 1].astype(int) - 1,
                         "low": rows[:, 2].astype(int) - 1, "rates": rows[:, 3:3 + n_temp]}
    return {"E": energy_cm * 1.4387769, "g": g,              # E in kelvin
            "up": rad[:, 1].astype(int) - 1, "low": rad[:, 2].astype(int) - 1,
            "A": rad[:, 3], "nu_ghz": rad[:, 4], "partners": partners}


def co_excitation_temperature(mol, t_cmb, n_h2, t_kin, n_levels=12):
    """T_ex of CO J=1-0 for an optically thin cloud bathed in the CMB.

    Rate equations for the lowest ``n_levels`` levels: spontaneous emission,
    absorption and stimulated emission of CMB photons (intensity B_nu(T_cmb)),
    and collisions with H2 at density n_h2 and temperature t_kin (thermal
    ortho/para ratio). With no collisions the levels follow the CMB exactly
    and T_ex = T_cmb; collisions push T_ex towards t_kin.
    """
    N = n_levels
    E, g = mol["E"][:N], mol["g"][:N]
    R = np.zeros((N, N))  # R[i, j]: rate from i to j
    for u, l, A, nu in zip(mol["up"], mol["low"], mol["A"], mol["nu_ghz"]):
        if u >= N:
            continue
        x = H_PLANCK * nu * 1e9 / (K_BOLTZ * t_cmb)
        occ = 1.0 / np.expm1(x)              # photon occupation number of the CMB
        R[u, l] += A * (1 + occ)              # spontaneous + stimulated emission
        R[l, u] += A * occ * g[u] / g[l]      # absorption
    opr = 9.0 * np.exp(-170.5 / t_kin)        # thermal ortho/para H2
    frac = {"oH2": opr / (1 + opr), "pH2": 1 / (1 + opr)}
    for key, part in mol["partners"].items():
        if key not in frac:
            continue
        k_ul = np.array([np.interp(t_kin, part["temps"], r) for r in part["rates"]])
        for u, l, k in zip(part["up"], part["low"], k_ul):
            if u >= N:
                continue
            c_ul = n_h2 * frac[key] * k
            R[u, l] += c_ul
            R[l, u] += c_ul * g[u] / g[l] * np.exp(-(E[u] - E[l]) / t_kin)
    M = R.T - np.diag(R.sum(axis=1))
    M[0, :] = 1.0                              # replace one equation by normalisation
    rhs = np.zeros(N)
    rhs[0] = 1.0
    pop = np.linalg.solve(M, rhs)
    dE = E[1] - E[0]
    return float(dE / np.log(pop[0] * g[1] / (pop[1] * g[0])))


def run():
    df = load_measurements()
    base = select(df, "baseline")
    out = {"n_baseline": int((base["kind"] != "upper").sum()),
           "n_upper_limits": int((base["kind"] == "upper").sum())}
    out["baseline"] = fit_beta(base)
    out["baseline_models"] = compare_models(base)
    out["t0_highz"] = fit_t0(base)
    by_method = {}
    for label, mask in (("SZ clusters", base["method"] == "tSZ"),
                        ("molecules and atoms", ~base["method"].isin(["tSZ", "CMB spectrum"]))):
        d = pd.concat([base[base["z"] == 0], base[mask]])
        by_method[label] = {**fit_beta(d), **compare_models(d)}
    out["by_method"] = by_method
    alt = pd.concat([base[~base["object"].str.startswith("PKS1830")],
                     select(df, "muller")])
    out["with_muller13"] = fit_beta(alt)
    # Replications.
    nt = select(df, "noterdaeme11")
    sz = nt["method"] == "tSZ"
    co = nt["method"].str.startswith("CO")
    carbon = nt["method"].str.startswith(("CI", "CII"))
    z0 = nt["z"] == 0
    out["noterdaeme11"] = {
        "SZ": fit_beta(nt[z0 | sz]), "SZ+C": fit_beta(nt[z0 | sz | carbon]),
        "SZ+CO": fit_beta(nt[z0 | sz | co]), "all": fit_beta(nt),
    }
    rch = select(df, "riechers22")
    out["riechers22_all_rows"] = fit_beta(rch)
    # One value per absorber: collision-corrected (Klimenko) where available.
    keep = ~rch["method"].str.contains("Maeder") & ~(
        rch["method"].isin(["CO Tex"]) & rch["object"].isin(
            rch.loc[rch["method"].str.contains("Klimenko") & rch["method"].str.startswith("CO"),
                    "object"]))
    out["riechers22_one_per_absorber"] = fit_beta(rch[keep])
    out["sz_null"] = observed_null(base)
    # CO thermometer physics.
    mol = read_lamda()
    dens = np.logspace(0, 4, 41)
    curves = {}
    for z in (0.0, 1.73, 2.69):
        t_cmb = T0_FIRAS[0] * (1 + z)
        curves[z] = np.array([co_excitation_temperature(mol, t_cmb, n, 100.0) for n in dens]) - t_cmb
    out["co_excitation"] = {"n_h2": dens, "t_kin": 100.0, "excess": curves}
    out["data"] = df
    return out
