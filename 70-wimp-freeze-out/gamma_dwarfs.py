"""Gamma rays from dwarf galaxies: our own combined limit from 14 years of Fermi-LAT data.

If WIMPs froze out by annihilating, they still annihilate today wherever dark matter is dense.
Each annihilation into b quarks or tau leptons ends in a shower of photons. The Milky Way's dwarf
spheroidal galaxies are the cleanest place to look: nearby, almost all dark matter, no gas and no
young stars to make gamma rays of their own. The expected photon flux from one dwarf is

    dPhi/dE = <sigma v> / (8 pi m^2) * dN/dE * J,       J = integral of rho^2 along the line of sight,

for a self-conjugate particle (Majorana); J is measured from the stars' velocity dispersions.

McDaniel et al. (2024) analysed 14.3 years of LAT data on 55 dwarfs and released, for each dwarf,
the likelihood of its photon counts as a function of the energy flux in each of 26 energy bins
from 0.5 GeV to 1 TeV (fermipy "SED" files). That is all a dark-matter search needs: for a mass
and a channel the predicted spectrum gives a flux in every bin, and the bin likelihoods add up.
We do that ourselves:

1. photons per annihilation from PPPC 4 DM ID (Cirelli et al. 2011), the tables McDaniel et al.
   used (CosmiXs as a check), integrated over each bin to an energy flux;
2. each dwarf's J-factor treated as a nuisance parameter with a log-normal prior (their Table 1)
   and profiled;
3. the profiles summed over a sample of dwarfs; the 95% upper limit where 2 Delta ln L = 2.71.

As a check we also redo their own simplified conversion (spectrum at the bin centre, as in their
release notebook) and compare it with the TS profiles they published, then sum their profiles.
"""

from __future__ import annotations

import ast
import json
import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from astropy.io import fits

import config as C

FERMI_DIR = C.RAW_DIR / "fermi"
MCD_DIR = FERMI_DIR / "mcdaniel2024"
CHANNEL_COLUMNS = {"pppc": {"bb": 13, "tau": 10}, "cosmixs": {"bb": 18, "tau": 10}}
LN10 = math.log(10.0)


# --- dwarf table and samples ----------------------------------------------------------------


@lru_cache(maxsize=None)
def dwarf_table() -> dict:
    import csv

    with open(FERMI_DIR / "mcdaniel2024_dSphs.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    return {r["name"]: {"log10J": float(r["log10_J"]), "sigmaJ": float(r["J_sigma"]),
                        "method": r["Method"], "distance_kpc": float(r["Distance"]),
                        "label": r["abbrv"], "ra": float(r["RA"]), "dec": float(r["Dec"])} for r in rows}


@lru_cache(maxsize=None)
def samples() -> dict:
    """The 'measured', 'benchmark' and 'inclusive' lists, read from the release's own notebook."""
    import warnings

    nb = json.loads((FERMI_DIR / "mcdaniel2024_basic_data_usage.ipynb").read_text())
    out = {}
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        with warnings.catch_warnings():          # their plot labels hold LaTeX like "$\\chi$"
            warnings.simplefilter("ignore", SyntaxWarning)
            tree = ast.parse("".join(cell["source"]).replace("%%time", ""))
        for node in tree.body:
            if (isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id in ("measured", "benchmark", "inclusive")
                    and isinstance(node.value, ast.List)):
                out[node.targets[0].id] = [ast.literal_eval(e) for e in node.value.elts]
    return out


@dataclass
class SED:
    name: str
    e_min: np.ndarray          # MeV
    e_max: np.ndarray
    e_ref: np.ndarray
    ref_eflux: np.ndarray      # MeV cm^-2 s^-1 (energy flux of the reference model in the bin)
    ref_dnde: np.ndarray       # MeV^-1 cm^-2 s^-1 at e_ref
    norm_scan: np.ndarray      # (bins, 20)
    dloglike: np.ndarray       # (bins, 20)
    ts: np.ndarray


@lru_cache(maxsize=None)
def load_sed(name: str) -> SED:
    with fits.open(MCD_DIR / "SEDs" / f"{name}_sed.fits") as h:
        d = h[1].data
        return SED(name, np.array(d["e_min"], float), np.array(d["e_max"], float), np.array(d["e_ref"], float),
                   np.array(d["ref_eflux"], float), np.array(d["ref_dnde"], float),
                   np.array(d["norm_scan"], float), np.array(d["dloglike_scan"], float), np.array(d["ts"], float))


def bin_loglike(sed: SED, norm: np.ndarray, beyond: str = "extrapolate") -> np.ndarray:
    """Sum over bins of the interpolated delta-log-likelihood at the given norms (bins, ...).

    Linear interpolation in norm. Beyond the scan's last point the last segment is extrapolated
    linearly ('extrapolate'), or the last value is kept ('clamp', what McDaniel et al.'s J-profile
    code does). The scans reach far past any limit, so this only matters far from the maximum."""
    total = np.zeros(norm.shape[1:])
    for j in range(norm.shape[0]):
        xs, ys = sed.norm_scan[j], sed.dloglike[j]
        nj = norm[j]
        v = np.interp(nj, xs, ys)
        hi = nj > xs[-1]
        if np.any(hi) and beyond == "extrapolate":
            slope = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
            v = np.where(hi, ys[-1] + slope * (nj - xs[-1]), v)
        total = total + v
    return total


# --- photons per annihilation -----------------------------------------------------------------


@lru_cache(maxsize=None)
def spectra_table(source: str = "pppc"):
    if source == "pppc":
        d = np.loadtxt(FERMI_DIR / "pppc4dmid_AtProduction_gammas.dat", skiprows=1)
    else:
        d = np.loadtxt(FERMI_DIR / "cosmixs_AtProduction-Gamma.dat", comments="#")
    masses = np.unique(d[:, 0])
    nx = int(np.sum(d[:, 0] == masses[0]))
    lx = d[:nx, 1]
    return masses, lx, d.reshape(len(masses), nx, d.shape[1])


def dn_dlog10x(m: float, channel: str, source: str = "pppc", mass_interp: str = "log"):
    """(log10 x grid, dN/dlog10x) for one mass: linear interpolation between tabulated masses in
    log m ('log', ours) or in m ('linear', as SciPy's interp2d in McDaniel et al.'s notebook)."""
    masses, lx, cube = spectra_table(source)
    col = CHANNEL_COLUMNS[source][channel]
    if not masses[0] <= m <= masses[-1]:
        raise ValueError(f"mass {m} outside the {source} table ({masses[0]}-{masses[-1]} GeV)")
    k = int(np.clip(np.searchsorted(masses, m) - 1, 0, len(masses) - 2))
    m0, m1 = masses[k], masses[k + 1]
    w = (math.log(m / m0) / math.log(m1 / m0)) if mass_interp == "log" else (m - m0) / (m1 - m0)
    return lx, (1 - w) * cube[k, :, col] + w * cube[k + 1, :, col]


def bin_energy_flux_per_annihilation(m: float, channel: str, e_min_mev, e_max_mev, source: str = "pppc"):
    """Integral of E dN/dE over each bin, in GeV per annihilation (E dN/dE dE = m x dN/dlog10x dlog10x)."""
    lx, dn = dn_dlog10x(m, channel, source)
    out = np.zeros(len(e_min_mev))
    for j, (lo, hi) in enumerate(zip(e_min_mev, e_max_mev)):
        a, b = math.log10(lo * 1e-3 / m), math.log10(hi * 1e-3 / m)
        if a >= 0:
            continue
        b = min(b, 0.0)
        g = np.linspace(a, b, 41)
        y = np.interp(g, lx, dn, left=0.0, right=0.0)
        out[j] = m * np.trapezoid(10 ** g * y, g)
    return out


def dnde_at(m: float, channel: str, e_mev: np.ndarray, source: str = "pppc", mass_interp: str = "linear"):
    """dN/dE (1/MeV) at energies e_mev, interpolated in log-log as McDaniel et al.'s notebook."""
    lx, dn = dn_dlog10x(m, channel, source, mass_interp)
    E_gev = 10 ** lx * m
    f = np.where(dn > 0, dn / (LN10 * E_gev), 1e-30)
    lf = np.log10(f)
    le = np.log10(e_mev * 1e-3)
    inside = (le >= np.log10(E_gev[0])) & (le <= np.log10(E_gev[-1]))
    vals = 10 ** np.interp(le, np.log10(E_gev), lf)
    return np.where(inside, vals, 0.0) * 1e-3


# --- likelihoods ---------------------------------------------------------------------------------

#: log10 sigma v grid (cm^3/s) for profiles and limits, and the nuisance grid in log10 J (in sigma).
LSV = np.arange(-28.5, -20.0001, 0.01)
NU = np.linspace(-5.0, 5.0, 201)


def dwarf_profile(name: str, m: float, channel: str, source: str = "pppc", method: str = "eflux",
                  prior: bool = True, lsv: np.ndarray = LSV) -> np.ndarray:
    """Profile log-likelihood (relative to no signal) on the log10 <sigma v> grid for one dwarf.

    method 'eflux': predicted energy flux in each bin over the reference energy flux (standard).
    method 'mcd': dN/dE at the bin's reference energy over ref_dnde (McDaniel et al.'s notebook).
    """
    sed = load_sed(name)
    dw = dwarf_table()[name]
    if method == "eflux":
        per_ann = bin_energy_flux_per_annihilation(m, channel, sed.e_min, sed.e_max, source) * 1e3  # MeV
        unit = per_ann / sed.ref_eflux          # norm per unit (sigma v J / 8 pi m^2)
        beyond, nu = "extrapolate", NU * dw["sigmaJ"]
    else:
        unit = dnde_at(m, channel, sed.e_ref, source) / sed.ref_dnde
        # their J profile: scale factor in [1e-3, 1e3] and the last scan value held beyond the scan
        beyond, nu = ("clamp" if prior else "extrapolate"), np.linspace(-3.0, 3.0, 601)
    l0 = float(bin_loglike(sed, np.zeros((len(unit), 1)))[0])
    if not prior:
        A = 10 ** (lsv + dw["log10J"]) / (8 * math.pi * m * m)
        return bin_loglike(sed, unit[:, None] * A[None, :], beyond) - l0
    A = 10 ** (lsv[:, None] + dw["log10J"] + nu[None, :]) / (8 * math.pi * m * m)
    ll = bin_loglike(sed, unit[:, None, None] * A[None, :, :], beyond) - l0
    ll -= 0.5 * (nu[None, :] / dw["sigmaJ"]) ** 2
    return ll.max(axis=1)


def upper_limit(lsv: np.ndarray, ts: np.ndarray, delta: float = C.DELTA_TS_95) -> dict:
    """95% one-sided limit: the largest <sigma v> at which TS is still within 'delta' of its maximum.

    For a profile that falls monotonically above the best fit this is simply where it crosses
    TS_max - 2.71. A stack can have a small bump instead (one dwarf's mild excess, e.g. Reticulum
    II's near 4e-25 cm^3/s at TeV masses, against the others' downward pull); taking the *last*
    crossing is then the conservative choice and keeps the limit continuous in mass."""
    i0 = int(np.argmax(ts))
    target = ts[i0] - delta
    inside = np.flatnonzero(ts[i0:] >= target)
    k = i0 + inside[-1]                      # last grid point still inside
    best = float(10 ** lsv[i0]) if i0 > 0 else 0.0
    if k + 1 >= len(ts):
        return {"ul": float("nan"), "best": best, "ts_max": float(ts[i0])}
    x = lsv[k] + (lsv[k + 1] - lsv[k]) * (ts[k] - target) / (ts[k] - ts[k + 1])
    return {"ul": float(10 ** x), "best": best, "ts_max": float(ts[i0])}


def combined_limits(sample: str, channel: str, masses, source: str = "pppc", method: str = "eflux",
                    prior: bool = True) -> dict:
    names = samples()[sample]
    out = {"mass": [], "ul": [], "ts_max": [], "best": [], "per_dwarf_ts_max": {}}
    for m in masses:
        tot = np.zeros_like(LSV)
        for n in names:
            ll = dwarf_profile(n, m, channel, source, method, prior)
            tot += 2 * ll
            out["per_dwarf_ts_max"].setdefault(n, []).append(float(2 * ll.max()))
        r = upper_limit(LSV, tot)
        out["mass"].append(float(m))
        out["ul"].append(r["ul"])
        out["ts_max"].append(r["ts_max"])
        out["best"].append(r["best"])
    return out


# --- McDaniel et al.'s own TS profiles -----------------------------------------------------------


def mcd_grid():
    m = np.logspace(math.log10(C.MCD_MASSES[0]), math.log10(C.MCD_MASSES[1]), C.MCD_MASSES[2])
    sv = np.logspace(math.log10(C.MCD_SIGMAV[0]), math.log10(C.MCD_SIGMAV[1]), C.MCD_SIGMAV[2])
    return m, sv


def mcd_ts(name: str, channel: str, prior: bool = True) -> np.ndarray:
    """Their published TS profile for one dwarf: array (40 masses, 60 cross sections)."""
    tag = "Jprior" if prior else "noprior"
    return np.load(MCD_DIR / "TS_profiles" / f"{name}_{tag}_{channel}.npy")


def mcd_stack_limits(sample: str, channel: str, prior: bool = True) -> dict:
    """Sum their TS profiles over a sample; limits by our upper_limit on a log-linear interpolation
    of their 60-point cross-section grid (their notebook interpolates linearly in sigma v)."""
    m, sv = mcd_grid()
    tot = sum(mcd_ts(n, channel, prior) for n in samples()[sample])
    ul, tsmax, best = [], [], []
    for i in range(len(m)):
        r = upper_limit(np.log10(sv), tot[i])
        ul.append(r["ul"])
        tsmax.append(r["ts_max"])
        best.append(r["best"])
    return {"mass": m.tolist(), "ul": ul, "ts_max": tsmax, "best": best}


def reproduce_mcd_profile(name: str, channel: str, masses_idx=(10, 20, 30), prior: bool = True) -> dict:
    """Recompute their TS profile from their SED with their conventions, at a few of their masses."""
    m, sv = mcd_grid()
    ref = mcd_ts(name, channel, prior)
    out = []
    for i in masses_idx:
        if m[i] < C.GAMMA_MASS_MIN:
            continue
        ll = dwarf_profile(name, float(m[i]), channel, "pppc", "mcd", prior, lsv=np.log10(sv))
        ts = 2 * ll
        near = ref[i] > ref[i].max() - 10.0            # where limits are decided
        out.append({"mass": float(m[i]), "max_abs_dTS_near_max": float(np.max(np.abs(ts - ref[i])[near])),
                    "max_abs_dTS_all": float(np.max(np.abs(ts - ref[i]))),
                    "ref_ts_max": float(ref[i].max()), "ours_ts_max": float(ts.max())})
    return {"dwarf": name, "channel": channel, "rows": out}


# --- cached per-dwarf profiles over a mass grid ------------------------------------------------

MASSES = np.logspace(math.log10(C.GAMMA_MASS_MIN), 4.0, 54)


def profile_cube(channel: str, source: str = "pppc", method: str = "eflux", force: bool = False) -> dict:
    """TS(<sigma v>) for every dwarf in any sample, every mass of MASSES; cached in data/derived."""
    path = C.DERIVED_DIR / f"dwarf_ts_{channel}_{source}_{method}.npz"
    names = sorted(set().union(*samples().values()))
    if path.exists() and not force:
        z = np.load(path)
        if list(z["names"]) == names and np.allclose(z["masses"], MASSES):
            return {"names": names, "masses": MASSES, "lsv": z["lsv"], "ts": z["ts"]}
    ts = np.zeros((len(names), len(MASSES), len(LSV)), dtype=np.float32)
    for i, n in enumerate(names):
        for j, m in enumerate(MASSES):
            ts[i, j] = 2 * dwarf_profile(n, float(m), channel, source, method, prior=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, names=np.array(names), masses=MASSES, lsv=LSV, ts=ts)
    return {"names": names, "masses": MASSES, "lsv": LSV, "ts": ts}


def stack_limits(cube: dict, sample: str) -> dict:
    idx = [cube["names"].index(n) for n in samples()[sample]]
    tot = cube["ts"][idx].astype(float).sum(axis=0)
    rows = [upper_limit(cube["lsv"], tot[j]) for j in range(len(cube["masses"]))]
    return {"mass": cube["masses"].tolist(), "ul": [r["ul"] for r in rows],
            "ts_max": [r["ts_max"] for r in rows], "best": [r["best"] for r in rows],
            "n_dwarfs": len(idx)}
