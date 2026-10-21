"""Chae's (2023, corrected in 2024a) acceleration-plane test, reimplemented.

For each wide binary the observed sky-projected separation s and velocity
v_p are "deprojected" to a 3D separation r and speed v by Monte Carlo: an
eccentricity is drawn from Hwang et al.'s (2022) individual estimate (a
split normal around their most likely value), the orbit is given a random
inclination (p ~ sin i), periastron angle and time, and then

    r = s / sqrt(cos^2 phi + cos^2 i sin^2 phi),
    v = v_p / sqrt(cos^2 psi + cos^2 i sin^2 psi),

with psi the direction of the orbital velocity in the orbit plane. Each
binary gives a point (g_N = G M / r^2, g = v^2 / r) in the "acceleration
plane", which in Newtonian dynamics scatters about g = g_N (2 - r/a), so
below the line for eccentric orbits. The test statistic is the median
offset Delta_perp = (log g - log g_N)/sqrt(2) in bins of
x0 = (log g_N + log g)/2, for the data and for a mock sample in which the
observed proper motions are replaced by Newtonian ones (same masses,
eccentricities, distances). Their difference delta_obs-newt, over 200 Monte
Carlo sets, measures the gravity boost: g_obs/g_pred = 10^(sqrt(2) delta).

Hidden close companions are given to a fraction f_multi of the binaries (40%
on star A only, 30% on B only, 30% on both): the star's light is split by a
magnitude difference with p(dM_G) ~ (dM_G/12)^-0.7 on 0-12, the masses come
from Chae's polynomial, the inner orbit has a log-flat semi-major axis from
0.01 AU to d AU (1" at the pair's distance) and its photocentre moves by
eta = M1 M2 (M1^2.5 - M2^2.5)/((M1 + M2)(M1^3.5 + M2^3.5)) times the inner
relative velocity, or not at all for inner periods under 3 yr.
f_multi is calibrated so that delta = 0 in the highest-acceleration bin.

Everything follows Chae's published code (Zenodo 10062232,
velprofile_wb_v3.py: Newton_simulation and Newton_2d_to_3d), vectorised.
"""

from __future__ import annotations

import numpy as np

from config import RAW_DIR

FAC = 1.0e-3 * (1 / (648000 / np.pi)) * 3.085677e13 / 3.15582e7   # pc mas/yr -> km/s (Chae's)
G_SI, MSUN, AU = 6.6743e-11, 1.9884e30, 1.496e11
X0_BINS = np.array([-11.5, -9.8, -8.5, -7.5])
ALPHA_L = 3.5


def chae_poly(which: str = "V"):
    """His M(M_G) and M_G(M) polynomial coefficients (log10 M = sum a_j M_G^j;
    M_G = sum b_j (log10 M)^j)."""
    tag = {"V": "V_Ic", "J": "J_Bp_Rp", "F": "F"}[which]
    d = RAW_DIR / "chae_zenodo" / "8416435"
    a = np.loadtxt(d / f"M_MagG_polyfit_{tag}.dat", skiprows=1, usecols=1) if (d / f"M_MagG_polyfit_{tag}.dat").exists() else None
    b = np.loadtxt(d / f"MagG_M_polyfit_{tag}.dat", skiprows=1, usecols=1)
    return a, b


def load_chae_base():
    import pandas as pd
    d = pd.read_csv(RAW_DIR / "chae_zenodo" / "8416435" / "gaia_dr3_MSMS_d200pc.csv")
    d.columns = [c.split("[")[0] for c in d.columns]
    return {c: d[c].to_numpy() for c in d.columns}


def chae_sample(base, *, dmax: float = 200.0, pmerr_max: float = 0.01, ext_corr: bool = True):
    """Chae (2023) main sample from his own base file: 26,615 pairs."""
    d = base
    MA = d["MagG1"] - (d["A_G_A"] if ext_corr else 0.0)
    MB = d["MagG2"] - (d["A_G_B"] if ext_corr else 0.0)
    a, b = 4 - 3.6 * 17 / 11.5 + 0.1, 3.6 / 11.5
    dec_min = -28.0 if ext_corr else -100.0
    m = (d["DEC1"] > dec_min) & (MA > 4) & (MA < 14) & (MB > 4) & (MB < 14)
    m &= (d["bp_rp1"] > a + b * MA) & (d["bp_rp2"] > a + b * MB)
    m &= d["R_chance"] < 0.01
    m &= (d["s"] > 0.2) & (d["s"] < 30)
    dd = np.abs(d["d1"] - d["d2"])
    m &= dd < np.sqrt((3 * d["d1_err"]) ** 2 + (3 * d["d2_err"]) ** 2 + (6 * d["s"] * 1e3 / 206265) ** 2)
    pmerr = np.maximum.reduce([np.abs(d[f"mu{w}{c}_err"] / d[f"mu{w}{c}"]) for w in (1, 2) for c in ("ra", "dec")])
    m &= pmerr < pmerr_max
    dM = (d["d1"] / d["d1_err"] ** 2 + d["d2"] / d["d2_err"] ** 2) / (1 / d["d1_err"] ** 2 + 1 / d["d2_err"] ** 2)
    m &= dM < dmax
    out = {k: v[m] for k, v in d.items()}
    out["MA"], out["MB"], out["dM"] = MA[m], MB[m], dM[m]
    return out


# --------------------------------------------------------------------------
# Kepler geometry
# --------------------------------------------------------------------------

def _true_anomaly(e, rng):
    M = 2 * np.pi * rng.random(e.shape)
    E = M + e * np.sin(M)
    for _ in range(50):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))


def _draw_e(e, e0, e1, rng):
    z = rng.standard_normal(e.shape)
    lo = np.maximum(e + z * (e - e0), 0.001)
    hi = np.minimum(e + z * (e1 - e), 0.999)
    return np.where(z < 0, lo, hi)


def _draw_e_stat(s_kau, rng):
    x = np.log10(s_kau * 1000)
    gam = np.clip(1 - 5.12340377 + 4.61235048 * x - 1.09835007 * x**2 + 0.08748435 * x**3, 1e-3, 2.5)
    return rng.random(s_kau.shape) ** (1.0 / gam)       # numpy power(a): pdf a x^(a-1)


def _angles(e, rng):
    inc = np.arccos(rng.random(e.shape))
    phi0 = 2 * np.pi * rng.random(e.shape)
    phi = phi0 + _true_anomaly(e, rng)
    psi = np.arctan(-(np.cos(phi) + e * np.cos(phi0)) / (np.sin(phi) + e * np.sin(phi0)))
    return inc, phi0, phi, psi


def deproject(s_kau, vp_kms, mtot, e, rng):
    """Chae's Newton_2d_to_3d: (log10 g_N, log10 g) in m/s^2."""
    inc, phi0, phi, psi = _angles(e, rng)
    r = s_kau / np.sqrt(np.cos(phi) ** 2 + np.cos(inc) ** 2 * np.sin(phi) ** 2)
    R = r * 1000 * AU
    gN = np.log10(G_SI * MSUN * mtot / R**2)
    v3 = vp_kms / np.sqrt(np.cos(psi) ** 2 + np.cos(inc) ** 2 * np.sin(psi) ** 2)
    g = np.log10((v3 * 1e3) ** 2 / R)
    return gN, g


# --------------------------------------------------------------------------
# Hidden companions (Chae 2023, Section 3.2)
# --------------------------------------------------------------------------

def companion_masses(MA, MB, f_multi, rng, *, gam_M: float = -0.7, poly: str = "V"):
    """Masses (A1, A2, B1, B2) with hidden companions for a random fraction
    f_multi (40% A only, 30% B only, 30% both)."""
    a, b = chae_poly(poly)
    mass = lambda MG: 10.0 ** np.polynomial.polynomial.polyval(MG, a)
    n = len(MA)
    has = rng.random(n) < f_multi
    which = rng.random(n)
    onA = has & ((which < 0.4) | (which >= 0.7))
    onB = has & (which >= 0.4)
    out = {}
    for name, MG, on in (("A", MA, onA), ("B", MB, onB)):
        dM = 12.0 * rng.random(n) ** (1.0 / (1.0 + gam_M))
        kappa = 1.0 / (1.0 + 10 ** (-0.4 * dM))
        Mh = np.where(on, -2.5 * np.log10(kappa) + MG, MG)
        Mc = -2.5 * np.log10(np.maximum(1 - kappa, 1e-20)) + MG
        out[name + "1"] = mass(Mh)
        out[name + "2"] = np.where(on, mass(Mc), 0.0)
    return out


def _inner_velocity(M1, M2, d_pc, rng, a_out_au):
    """Photocentre velocity (km/s, x/y) of a hidden inner pair (Chae's code)."""
    n = len(M1)
    lo = np.log10(0.01)
    hi = np.log10(np.minimum(d_pc, a_out_au * 1000 * 0.3))
    a_in = 10 ** (lo + (hi - lo) * rng.random(n))
    gam = np.clip(1 - 1.26411333 + 0.85366141 * np.log10(a_in), 0.001, 2.3)
    e = rng.random(n) ** (1.0 / gam)
    inc, phi0, phi, psi = _angles(e, rng)
    r_in = a_in * (1 - e**2) / (1 + e * np.cos(phi - phi0))
    P_in = np.sqrt(a_in**3 / (M1 + M2))
    eta = np.where(P_in < 3, 0.0, M1 * M2 * (M1 ** (ALPHA_L - 1) - M2 ** (ALPHA_L - 1))
                   / ((M1 + M2) * (M1**ALPHA_L + M2**ALPHA_L)))
    v = eta * 1e-3 * np.sqrt((G_SI * MSUN * (M1 + M2) / AU) * np.maximum(2 / r_in - 1 / a_in, 0))
    return v * np.cos(psi), v * np.cos(inc) * np.sin(psi)


def newton_mock_vp(d, m, e, rng, gboost: float = 1.0, noise: bool = False, persp=None):
    """Chae's Newton_simulation: mock sky-projected velocities (km/s).

    ``noise``: add Gaussian proper-motion errors (each star's catalogue
    errors) to the mock, as the real proper motions carry them (Chae's code
    compares noise-free mock velocities with the measured, noisy ones).
    ``persp``: (dmu_ra, dmu_dec) perspective term per pair to add likewise."""
    s = d["s"]
    MA, MB = m["A1"] + m["A2"], m["B1"] + m["B2"]
    inc, phi0, phi, psi = _angles(e, rng)
    r = 1e3 * s / np.sqrt(np.cos(phi) ** 2 + np.cos(inc) ** 2 * np.sin(phi) ** 2)   # AU
    a = r * (1 + e * np.cos(phi - phi0)) / (1 - e**2)
    vc = 1e-3 * np.sqrt((G_SI * gboost * MSUN * (MA + MB) / AU) * np.maximum(2 / r - 1 / a, 0))
    dr_los = r * np.sin(inc) * np.sin(phi) / 206265
    w1, w2 = 1 / d["d1_err"] ** 2, 1 / d["d2_err"] ** 2
    dM0 = (d["d1"] * w1 + d["d2"] * w2) / (w1 + w2)
    dA = dM0 + MB / (MA + MB) * dr_los
    dB = dM0 - MA / (MA + MB) * dr_los
    dM = (dA * w1 + dB * w2) / (w1 + w2)
    dvx, dvy = vc * np.cos(psi), vc * np.cos(inc) * np.sin(psi)
    muxM = 0.5 * (d["mu1ra"] + d["mu2ra"])
    muyM = 0.5 * (d["mu1dec"] + d["mu2dec"])
    muxA = muxM + MB / (MA + MB) * dvx / (FAC * dA)
    muxB = muxM - MA / (MA + MB) * dvx / (FAC * dB)
    muyA = muyM + MB / (MA + MB) * dvy / (FAC * dA)
    muyB = muyM - MA / (MA + MB) * dvy / (FAC * dB)
    for X, dX in (("A", dA), ("B", dB)):
        has = m[X + "2"] > 0
        if np.any(has):
            vx, vy = _inner_velocity(m[X + "1"][has], m[X + "2"][has], dM[has], rng, a[has])
            if X == "A":
                muxA[has] += vx / (FAC * dX[has])
                muyA[has] += vy / (FAC * dX[has])
            else:
                muxB[has] += vx / (FAC * dX[has])
                muyB[has] += vy / (FAC * dX[has])
    if noise:
        muxA = muxA + d["mu1ra_err"] * rng.standard_normal(len(s))
        muyA = muyA + d["mu1dec_err"] * rng.standard_normal(len(s))
        muxB = muxB + d["mu2ra_err"] * rng.standard_normal(len(s))
        muyB = muyB + d["mu2dec_err"] * rng.standard_normal(len(s))
    if persp is not None:      # spurious relative PM (star B minus A) from the systemic motion
        muxB = muxB + persp[0]
        muyB = muyB + persp[1]
    return FAC * np.hypot(muxA - muxB, muyA - muyB) * dM


def observed_vp(d):
    w1, w2 = 1 / d["d1_err"] ** 2, 1 / d["d2_err"] ** 2
    dM = (d["d1"] * w1 + d["d2"] * w2) / (w1 + w2)
    return FAC * np.hypot(d["mu1ra"] - d["mu2ra"], d["mu1dec"] - d["mu2dec"]) * dM


def binned_medians(gN, g, bins=X0_BINS):
    x0 = 0.5 * (gN + g)
    dperp = (g - gN) / np.sqrt(2)
    idx = np.digitize(x0, bins) - 1
    return np.array([np.median(dperp[idx == k]) if np.any(idx == k) else np.nan for k in range(len(bins) - 1)])


def run(d, f_multi: float, *, n_mc: int = 200, seed: int = 0, echoice: str = "i", poly: str = "V",
        gboost: float = 1.0, mock_noise: bool = False, mock_persp: bool = False):
    """delta_obs-newt in the three x0 bins for n_mc Monte Carlo sets.
    ``mock_noise``/``mock_persp`` give the Newtonian mock the measurement
    errors and perspective effects that the observed velocities contain."""
    rng = np.random.default_rng(seed)
    persp = None
    if mock_persp:
        import injection as INJ
        persp = INJ.perspective_dmu(d, rng=np.random.default_rng(seed + 99))
    vp_obs = observed_vp(d)
    out = np.empty((n_mc, len(X0_BINS) - 1))
    med_obs = np.empty_like(out)
    med_new = np.empty_like(out)
    for k in range(n_mc):
        m = companion_masses(d["MA"], d["MB"], f_multi, rng, poly=poly)
        mtot = m["A1"] + m["A2"] + m["B1"] + m["B2"]
        draw = (lambda: _draw_e(d["e"], d["e0"], d["e1"], rng)) if echoice == "i" else \
               (lambda: _draw_e_stat(d["s"], rng))
        gN, g = deproject(d["s"], vp_obs, mtot, draw(), rng)
        vp_mock = newton_mock_vp(d, m, draw(), rng, gboost=gboost, noise=mock_noise, persp=persp)
        gN2, g2 = deproject(d["s"], vp_mock, mtot, draw(), rng)
        med_obs[k] = binned_medians(gN, g)
        med_new[k] = binned_medians(gN2, g2)
        out[k] = med_obs[k] - med_new[k]
    return out, med_obs, med_new


def calibrate_fmulti(d, *, grid=None, n_mc: int = 60, seed: int = 0, **kw):
    """f_multi that zeroes delta in the highest-acceleration bin (x0 ~ -8)."""
    grid = np.linspace(0.0, 1.0, 11) if grid is None else grid
    top = []
    for f in grid:
        dl, _, _ = run(d, f, n_mc=n_mc, seed=seed, **kw)
        top.append(np.median(dl[:, -1]))
    top = np.array(top)
    # delta_top decreases with f_multi (more hidden mass -> higher Newtonian prediction)
    if np.all(top > 0):
        return grid[-1], grid, top
    if np.all(top < 0):
        return grid[0], grid, top
    j = np.flatnonzero(np.diff(np.sign(top)) != 0)[0]
    f = grid[j] + (0 - top[j]) * (grid[j + 1] - grid[j]) / (top[j + 1] - top[j])
    return float(f), grid, top
