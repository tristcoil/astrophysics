"""The Galaxy's double white dwarfs as LISA will hear them: a model built from observations.

Ingredients, each measured rather than simulated:

1. How many white dwarfs there are near the Sun: the volume-complete 40 pc sample
   of O'Brien et al. (2024), every one spectroscopically confirmed (B10 will study
   these stars in their own right; here we only count them and use their masses).
2. How the Galaxy's stars are spread out: the McMillan (2017) mass model (thin and
   thick exponential discs and a bulge). White dwarfs are assumed to trace stellar
   mass, so the local ratio (white dwarfs per solar mass of stars) scales the 40 pc
   count to the whole Galaxy.
3. How many white dwarfs have a close white-dwarf companion, and at what separations:
   Maoz, Hallakoun & Badenes (2018), from the radial-velocity variability of ~4,000
   SDSS and 439 SPY white dwarfs: 9.5% of white dwarfs are in pairs closer than 4 au,
   born with separations dN/da ~ a^alpha, alpha = -1.30, then shrunk by gravitational
   waves over 10 Gyr of constant star formation, which gives the present-day
   distribution (Maoz et al. 2012):
       N(x) ~ x^(4+alpha) [1 - (1 + x^-4)^((alpha+1)/4)],   x = a / (K t0)^(1/4),
       K = (256/5) G^3 m1 m2 (m1 + m2) / c^5.
4. Masses: the more massive star drawn from the 40 pc white-dwarf masses, its
   companion uniformly between 0.15 Msun (the lightest known He-core white dwarfs)
   and the primary (Korol et al. 2022's prescription; the RV surveys cannot pin the
   mass ratio).

This is the recipe of Korol et al. (2022), with the 40 pc sample in place of the
20 pc density and the SDSS mass function, and the full McMillan Galaxy in place of a
single disc. ``variant="korol22"`` reproduces their setup as a check.

LISA's view: every binary with f > 0.1 mHz gets a GW amplitude from its distance and
chirp mass and a random orientation; the confusion foreground is found by the usual
iteration (e.g. Karnesis et al. 2021): estimate the noise from all sources not yet
resolved, call a source resolved if its 4-year SNR exceeds 7 against instrument +
foreground noise, remove it, repeat until nothing changes. The SNRs are LISA's
sky-averaged ones (gw.py), so this is a sky-averaged model of the foreground.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C
import fetch_data as F
import gw

AU = 1.495978707e11
A_MIN = 2e7  # m: two white dwarfs touching (Maoz et al. 2012)
A_MAX = 4 * AU
T0_YR = 10e9
F_MIN = 1e-4  # Hz, LISA's nominal low-frequency end
F_MAX = 0.1

MAOZ18 = {"f_bin": (0.095, 0.020), "alpha": (-1.30, 0.15), "R_merge_per_wd_yr": (9.7e-12, 1.1e-12)}
#: McMillan (2017) best-fitting stellar components (as implemented in galpy's McMillan17).
MCMILLAN17 = {
    "R0_kpc": 8.21, "z_sun_kpc": 0.0208,
    "thin": {"Sigma0": 896.0, "Rd": 2.5, "zd": 0.3},  # Msun/pc^2, kpc, kpc
    "thick": {"Sigma0": 183.0, "Rd": 3.02, "zd": 0.9},
    "bulge": {"rho0": 98.4, "r0": 0.075, "rcut": 2.1, "alpha": 1.8, "q": 0.5},  # Msun/pc^3, kpc
}
#: Korol et al. (2022): one exponential disc with a sech^2 vertical profile, cut at 20 kpc.
KOROL22 = {"R0_kpc": 8.1, "z_sun_kpc": 0.03, "Rd": 2.5, "zd": 0.3, "Rmax": 20.0,
           "rho_wd_sun": 4.49e-3, "N_wd_paper": 2.7e9, "N_band_paper": 26e6, "N_detected_paper": 60e3}
#: Kepler et al. (2015) SDSS white-dwarf mass function (as used by Korol et al. 2022).
KEPLER15 = {"mu": [0.65, 0.57, 0.81], "sigma": [0.044, 0.097, 0.187], "w": [0.81, 0.14, 0.05]}


# ---------------------------------------------------------------------------
# Local white dwarfs
# ---------------------------------------------------------------------------


def read_40pc() -> pd.DataFrame:
    df = pd.read_csv(F.OBRIEN_40PC, sep="\t", comment="#")
    df = df.iloc[2:].reset_index(drop=True)
    for c in ("Plx", "e_Plx", "Mass", "e_Mass", "Massc", "Teff"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def local_density(df: pd.DataFrame | None = None) -> dict:
    """White dwarfs per pc^3 within 40 pc (Poisson error), and their mass sample."""
    df = read_40pc() if df is None else df
    inside = df[df["Plx"] >= 25.0]
    n = len(inside)
    vol = 4 / 3 * np.pi * 40.0**3
    masses = inside["Mass"].dropna().to_numpy()
    return {"n_wd": n, "volume_pc3": vol, "density": n / vol, "density_err": np.sqrt(n) / vol,
            "n_with_mass": int(masses.size), "masses": masses,
            "mass_median": float(np.median(masses)), "frac_below_045": float(np.mean(masses < 0.45))}


# ---------------------------------------------------------------------------
# Galaxy models
# ---------------------------------------------------------------------------


def _bulge_rho(m, b):
    return b["rho0"] * (1 + m / b["r0"]) ** (-b["alpha"]) * np.exp(-((m / b["rcut"]) ** 2))


def mcmillan_masses():
    g = MCMILLAN17
    out = {}
    for comp in ("thin", "thick"):
        c = g[comp]
        out[comp] = 2 * np.pi * c["Sigma0"] * (c["Rd"] * 1e3) ** 2  # Msun
    b = g["bulge"]
    m = np.linspace(0, 20 * b["rcut"], 200001)
    out["bulge"] = float(np.trapezoid(4 * np.pi * b["q"] * _bulge_rho(m, b) * m**2, m) * 1e9)  # kpc^3 -> pc^3
    return out


def mcmillan_local_stellar_density():
    g = MCMILLAN17
    R, z = g["R0_kpc"], g["z_sun_kpc"]
    rho = 0.0
    for comp in ("thin", "thick"):
        c = g[comp]
        rho += c["Sigma0"] / (2 * c["zd"] * 1e3) * np.exp(-R / c["Rd"] - abs(z) / c["zd"])
    b = g["bulge"]
    rho += _bulge_rho(np.hypot(R, z / b["q"]), b)
    return rho  # Msun/pc^3


def sample_positions(n, rng, variant):
    """Galactocentric (x, y, z) in kpc for n white dwarfs."""
    if variant == "korol22":
        k = KOROL22
        R = rng.gamma(2.0, k["Rd"], size=int(n * 1.2))
        R = R[R <= k["Rmax"]][:n]
        while R.size < n:
            extra = rng.gamma(2.0, k["Rd"], size=n)
            R = np.concatenate([R, extra[extra <= k["Rmax"]]])[:n]
        z = k["zd"] * np.arctanh(rng.uniform(-1, 1, size=n)) + 0.0  # sech^2 profile
        phi = rng.uniform(0, 2 * np.pi, size=n)
        return R * np.cos(phi), R * np.sin(phi), z
    masses = mcmillan_masses()
    tot = sum(masses.values())
    comp = rng.choice(3, size=n, p=[masses["thin"] / tot, masses["thick"] / tot, masses["bulge"] / tot])
    x, y, z = np.empty(n), np.empty(n), np.empty(n)
    for i, name in enumerate(("thin", "thick")):
        sel = comp == i
        c = MCMILLAN17[name]
        R = rng.gamma(2.0, c["Rd"], size=sel.sum())
        phi = rng.uniform(0, 2 * np.pi, size=sel.sum())
        x[sel], y[sel] = R * np.cos(phi), R * np.sin(phi)
        z[sel] = rng.laplace(0.0, c["zd"], size=sel.sum())
    sel = comp == 2
    b = MCMILLAN17["bulge"]
    mg = np.linspace(1e-5, 8 * b["rcut"], 40001)
    cdf = np.cumsum(_bulge_rho(mg, b) * mg**2)
    cdf /= cdf[-1]
    m = np.interp(rng.uniform(size=sel.sum()), cdf, mg)
    u = rng.uniform(-1, 1, size=sel.sum())
    phi = rng.uniform(0, 2 * np.pi, size=sel.sum())
    s = np.sqrt(1 - u**2)
    x[sel], y[sel], z[sel] = m * s * np.cos(phi), m * s * np.sin(phi), m * u * b["q"]
    return x, y, z


def total_white_dwarfs(rho_wd_sun, variant):
    if variant == "korol22":
        k = KOROL22
        # integral of exp(-(R - R0)/Rd) sech^2(z/zd) over the disc (R < Rmax), z_sun offset ignored
        Rg = np.linspace(0, k["Rmax"], 20001)
        surf = np.trapezoid(np.exp(-(Rg - k["R0_kpc"]) / k["Rd"]) * 2 * np.pi * Rg, Rg) * 1e6  # pc^2
        return rho_wd_sun * surf * 2 * k["zd"] * 1e3
    masses = mcmillan_masses()
    return rho_wd_sun / mcmillan_local_stellar_density() * sum(masses.values())


def sun_position(variant):
    g = KOROL22 if variant == "korol22" else MCMILLAN17
    return np.array([g["R0_kpc"], 0.0, g["z_sun_kpc"]])


# ---------------------------------------------------------------------------
# Separations (Maoz et al. 2012, 2018)
# ---------------------------------------------------------------------------


def n_of_x(x, alpha):
    x = np.asarray(x, float)
    if abs(alpha + 1) < 1e-9:
        return x**3 * np.log1p(x**-4.0)
    s = (alpha + 1) / 4
    return np.abs(x ** (4 + alpha) * (1 - (1 + x**-4.0) ** s))


class SeparationTable:
    def __init__(self, alpha):
        self.alpha = alpha
        self.lx = np.linspace(-6, 6, 120001)
        x = 10**self.lx
        dens = n_of_x(x, alpha) * x * np.log(10)  # per unit log10 x
        cdf = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(self.lx))])
        self.cdf = cdf

    def F(self, x):
        return np.interp(np.log10(x), self.lx, self.cdf)

    def inv(self, u):
        return 10 ** np.interp(u, self.cdf, self.lx)

    def density(self, x):
        return n_of_x(x, self.alpha)


def k_coefficient(m1, m2):
    return 256 / 5 * C.GM_SUN**3 * m1 * m2 * (m1 + m2) / C.C_LIGHT**5  # m^4/s


def sample_masses(n, rng, primary="40pc", masses_40pc=None):
    if primary == "kepler15":
        k = KEPLER15
        comp = rng.choice(3, size=n, p=np.array(k["w"]) / np.sum(k["w"]))
        m1 = rng.normal(np.array(k["mu"])[comp], np.array(k["sigma"])[comp])
    else:
        m1 = rng.choice(masses_40pc, size=n)
    m1 = np.clip(m1, 0.16, 1.38)
    m2 = rng.uniform(0.15, m1)
    elm = m1 < 0.25  # Korol et al.: very light primaries get a 0.2-1.2 Msun partner
    m2[elm] = rng.uniform(0.2, 1.2, size=elm.sum())
    return np.maximum(m1, m2), np.minimum(m1, m2)


# ---------------------------------------------------------------------------
# The population in the LISA band
# ---------------------------------------------------------------------------


def simulate(variant="mcmillan17", primary="40pc", alpha=MAOZ18["alpha"][0], f_bin=MAOZ18["f_bin"][0],
             rho_wd_sun=None, t0_yr=T0_YR, seed=2026, n_mass_draws=2_000_000, chunk=2_000_000,
             keep_sample=200_000):
    """Every DWD with f > F_MIN as (f, A^2 * inclination factor), plus summary numbers."""
    rng = np.random.default_rng(seed)
    loc = local_density()
    if rho_wd_sun is None:
        rho_wd_sun = KOROL22["rho_wd_sun"] if variant == "korol22" else loc["density"]
    n_wd = total_white_dwarfs(rho_wd_sun, variant)
    table = SeparationTable(alpha)
    t0 = t0_yr * C.YEAR
    # fraction of DWDs (a < 4 au) that are in the LISA band, averaged over masses
    m1, m2 = sample_masses(n_mass_draws, rng, primary, loc["masses"])
    scale = (k_coefficient(m1, m2) * t0) ** 0.25
    a_band = gw.semi_major_axis(2.0 / F_MIN, m1 + m2)
    F_min, F_band, F_max = table.F(A_MIN / scale), table.F(a_band / scale), table.F(A_MAX / scale)
    p_band = np.clip((F_band - F_min) / (F_max - F_min), 0, 1)
    frac_band = float(np.mean(p_band))
    n_dwd = n_wd * f_bin
    n_band = int(round(n_dwd * frac_band))
    # merger rate per WD in this model. Maoz et al. count the present-day pairs whose
    # merger time a^4/K is shorter than some interval, per interval; for 10 Myr that is the
    # number to compare with their 9.7e-12/yr. The instantaneous flux through contact is
    # ~3x higher: the power law of birth separations runs down to contact, so many pairs
    # are born almost touching and merge within years (a feature of the model, not data).
    x_min = A_MIN / scale
    dn_da = table.density(x_min) / scale / (F_max - F_min)
    rate_per_dwd = dn_da * k_coefficient(m1, m2) / (4 * A_MIN**3) * C.YEAR
    r_merge_instant = float(np.mean(rate_per_dwd)) * f_bin
    window_yr = 1e7
    x_win = np.maximum((window_yr / t0_yr) ** 0.25, x_min)
    r_merge = float(np.mean((table.F(x_win) - F_min) / (F_max - F_min))) * f_bin / window_yr
    # draw the in-band systems in chunks
    sun = sun_position(variant)
    w = p_band / p_band.sum()
    f_all, pw_all = [], []
    keep = {k: [] for k in ("f", "mc", "d_kpc", "m1", "m2", "A")}
    done = 0
    while done < n_band:
        n = min(chunk, n_band - done)
        idx = rng.choice(m1.size, size=n, p=w)
        mm1, mm2, sc = m1[idx], m2[idx], scale[idx]
        u = rng.uniform(F_min[idx], F_band[idx])
        a = table.inv(u) * sc
        P = 2 * np.pi * np.sqrt(a**3 / (C.GM_SUN * (mm1 + mm2)))
        f = 2.0 / P
        x, y, z = sample_positions(n, rng, variant)
        d = np.sqrt((x - sun[0]) ** 2 + (y - sun[1]) ** 2 + (z - sun[2]) ** 2)
        d = np.maximum(d, 1e-3)  # 1 pc floor (a handful of draws)
        mc = gw.chirp_mass(mm1, mm2)
        A = gw.amplitude(f, mc, d * C.KPC)
        fac = gw.inclination_factor(np.arccos(rng.uniform(-1, 1, size=n)))
        f_all.append(f.astype(np.float64))
        pw_all.append(A**2 * fac)  # float64: A^2 ~ 1e-46 underflows float32
        if sum(len(v) for v in keep["f"]) < keep_sample:
            k = slice(0, min(n, keep_sample))
            for name, arr in (("f", f), ("mc", mc), ("d_kpc", d), ("m1", mm1), ("m2", mm2), ("A", A)):
                keep[name].append(arr[k])
        done += n
    summary = {"variant": variant, "primary_masses": primary, "alpha": alpha, "f_bin": f_bin,
               "t0_yr": t0_yr, "rho_wd_sun_pc3": rho_wd_sun, "n_wd_galaxy": n_wd, "n_dwd_4au": n_dwd,
               "frac_in_band": frac_band, "n_in_band": n_band, "merger_rate_per_wd_yr": r_merge,
               "merger_rate_window_yr": window_yr, "merger_rate_per_wd_instantaneous_yr": r_merge_instant}
    sample = {k: np.concatenate(v)[:keep_sample] for k, v in keep.items()}
    return np.concatenate(f_all), np.concatenate(pw_all), summary, sample


def _block_median(values, block):
    n = values.size // block * block
    return np.median(values[:n].reshape(-1, block), axis=1)


def resolve(f, power, t_obs_yr=4.0, snr_threshold=7.0, block=1024, max_iter=40):
    """Iterative subtraction: which sources stand out above instrument + foreground noise.

    ``power`` is A^2 x inclination factor. A source's mean-square signal in the
    sky-averaged detector is 2 A^2 x factor, all of it in one frequency bin of width
    1/T, so the unresolved sources make a power spectrum P_k = T sum(2 A^2 factor)
    over the sources in bin k. As in the LISA pipelines (Karnesis et al. 2021; Korol
    et al. 2022), the foreground level is the running *median* of that spectrum
    (here: medians of blocks of ``block`` bins, ~8 microhertz for 4 years), which a few
    bright isolated sources cannot inflate -- above ~2 mHz most bins hold no source
    at all. SNR^2 = 4 T A^2 factor / (S_instrument + S_foreground).
    """
    T = t_obs_yr * C.YEAR
    df = 1.0 / T
    k = np.floor((f - F_MIN) / df).astype(np.int64)
    n_bins = int(np.ceil((F_MAX - F_MIN) / df)) + 1
    k = np.clip(k, 0, n_bins - 1)
    centres = F_MIN + (np.arange(n_bins // block) + 0.5) * block * df
    s_inst_c = gw.sn_instrument(centres)
    resolved = np.zeros(f.size, bool)
    history = []
    w2 = 2.0 * T * power.astype(np.float64)

    def foreground(mask):
        spec = np.bincount(k[mask], weights=w2[mask], minlength=n_bins)
        return _block_median(spec, block)

    for _ in range(max_iter):
        s_c = foreground(~resolved)
        snr2 = 4 * T * power / np.interp(f, centres, s_inst_c + s_c)
        new = snr2 > snr_threshold**2
        history.append(int(new.sum()))
        if np.array_equal(new, resolved):
            break
        resolved = new
    s_c = foreground(~resolved)
    spec_all = np.bincount(k, weights=w2, minlength=n_bins)
    s_all_mean = spec_all[: centres.size * block].reshape(-1, block).mean(axis=1)
    counts = np.bincount(k, minlength=n_bins)[: centres.size * block].reshape(-1, block).mean(axis=1)
    above = centres[s_c > s_inst_c]
    return {"f": centres, "s_conf": s_c, "s_all_sources_mean": s_all_mean, "s_inst": s_inst_c,
            "n_resolved": int(resolved.sum()), "iterations": history, "resolved": resolved,
            "f_conf_edge": float(above.max()) if above.size else np.nan, "n_per_bin_1overT": counts}
