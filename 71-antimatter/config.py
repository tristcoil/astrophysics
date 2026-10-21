"""Project configuration: paths, data sources, physical constants, the
cosmology we adopt, and the published values our results are checked against.

Published numbers are kept in one place (``PUBLISHED``) and are used only for
comparison, never silently as inputs. Where a published number *is* an input
(a measured cross-section, a lattice result, an experiment's upper limit),
it is listed under ``INPUTS`` with its source, and the README says so.
"""

from __future__ import annotations

import math
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CRDB_DIR = RAW_DIR / "crdb"
ARXIV_DIR = RAW_DIR / "arxiv"
EXTERNAL_DIR = DATA_DIR / "external"
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, CRDB_DIR, ARXIV_DIR, EXTERNAL_DIR, DERIVED_DIR, FIGURE_DIR,
           RAW_DIR / "vizier", RAW_DIR / "pdg"):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astro-atlas-P5/1.0 (python-urllib; research use)"

# ---------------------------------------------------------------------------
# Earlier atlas projects this one builds on
# ---------------------------------------------------------------------------
#
# D3 measured the baryon density from the Planck temperature spectrum; we read
# its results file. P2 (being built in parallel) downloaded Saikawa & Shirai's
# g*(T) table; we hard-link that file after a checksum check, or download it.

D3_DIR = ROOT_DIR / "09-planck-cmb"
D3_RESULTS = D3_DIR / "data" / "results.json"
#: D3's TT-only CAMB fit, as quoted in its README (used only if the results
#: file is missing).
D3_OMBH2_FALLBACK = (0.02221, 0.00019)

P2_DIR = ROOT_DIR / "68-thermal-history"

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: Cosmic-Ray DataBase (Maurin, Melot & Taillet 2014; Maurin et al. 2020, 2023).
#: Its REST interface returns CSV; combo_level=1 lets it form exact ratios from
#: fluxes published in the same paper (AMS-02's 2025 p-bar/p).
CRDB_URL = "http://lpsc.in2p3.fr/crdb/rest.php"
CRDB_QUERIES = {
    # file stem: (numerator, denominator, combo_level, energy axis). Data are
    # requested on their native axis (no conversion by CRDB): AMS-02 publishes
    # versus rigidity R, PAMELA and BESS versus kinetic energy per nucleon.
    "pbar_over_p": ("1H-bar", "1H", 1, "R"),
    "pbar_over_p_ekn": ("1H-bar", "1H", 1, "EKN"),
    "pbar_flux": ("1H-bar", "", 0, "R"),
    "p_flux": ("1H", "", 0, "R"),            # AMS-02 2025 (to 100 GV)
    "h_flux": ("H", "", 0, "R"),             # AMS-02 2015, 2021 (to 1.8 TV): filed as 'H'
    "he_flux": ("He", "", 0, "R"),
    "c_flux": ("C", "", 0, "R"),
    "o_flux": ("O", "", 0, "R"),
    "ne_flux": ("Ne", "", 0, "R"),
    "mg_flux": ("Mg", "", 0, "R"),
    "si_flux": ("Si", "", 0, "R"),
    "n_flux": ("N", "", 0, "R"),
    "b_over_c": ("B", "C", 0, "R"),
    "hebar_over_he": ("He-bar", "He", 0, "R"),
}

#: Saikawa & Shirai 2018 (JCAP 05, 035): g*rho(T), g*s(T) of the Standard Model
#: with lattice QCD (also used by P2).
SS2018_URL = "https://member.ipmu.jp/satoshi.shirai/standardmodel2018.dat"
SS2018_FILE = RAW_DIR / "ss2018_standardmodel.dat"
SS2018_SHA256 = "268a51b869198b13aef052185323e0b0015f78f4a0902f19a79e71d1017dfaad"
SIBLING_COPIES = {SS2018_FILE: P2_DIR / "data" / "raw" / "ss2018_standardmodel.dat"}

#: Fermi-LAT isotropic gamma-ray background (Ackermann et al. 2015, ApJ 799,
#: 86): the journal's machine-readable Table 3 on VizieR, 26 bands x Galactic
#: foreground models A, B, C; band-integrated intensities (ph cm^-2 s^-1 sr^-1).
VIZIER_DIR = RAW_DIR / "vizier"
VIZIER_FILES = {
    "J_ApJ_799_86_table3.dat": "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/799/86/table3.dat",
    "J_ApJ_799_86_ReadMe": "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/799/86/ReadMe",
}

#: PRyMordial (Burns, Tait & Valli 2024, EPJC 84, 86), the BBN code we run for
#: deuterium; GPL-3.0, not on PyPI. Pinned to the 2026-05-08 commit that
#: replaced the NACRE-II d(p,gamma)3He table by a thermal average of the LUNA
#: (Mossa et al. 2020) fit. Unpacked under data/external/ (not committed).
PRYMORDIAL_COMMIT = "725d8a8db3ad5ea2630580d825c9d0d69ed76533"
PRYMORDIAL_URL = f"https://codeload.github.com/vallima/PRyMordial/zip/{PRYMORDIAL_COMMIT}"
PRYMORDIAL_DIR = EXTERNAL_DIR / f"PRyMordial-{PRYMORDIAL_COMMIT}"

#: PDG compilations of measured antiproton-proton total and elastic cross
#: sections (rpp2022 files, unchanged in the 2024-2026 editions).
PDG_DIR = RAW_DIR / "pdg"
PDG_FILES = {
    "rpp2022-pbarp_total.dat": "https://pdg.lbl.gov/2026/hadronic-xsections/rpp2022-pbarp_total.dat",
    "rpp2022-pbarp_elastic.dat": "https://pdg.lbl.gov/2026/hadronic-xsections/rpp2022-pbarp_elastic.dat",
}

#: arXiv sources: we keep only the files we read (figure PostScript for
#: curve extraction, TeX for tables). id -> {local name: name in tarball}.
ARXIV_SOURCES = {
    # Cohen, De Rujula & Glashow 1998 (ApJ 495, 539): the annihilation rate
    # J(y), the plasma temperature and their predicted gamma-ray spectra, as
    # vector PostScript (XMgr and Mathematica) -- our answer key for the
    # domain calculation.
    "astro-ph/9707087": {"cdg98_J.eps": "curr.eps", "cdg98_T.eps": "temp.eps",
                         "cdg98_flux.eps": "flux.eps", "cdg98_dist.eps": "dist.eps",
                         "cdg98.tex": "final.tex"},
}

# ---------------------------------------------------------------------------
# Physical constants (SI unless stated; CODATA via scipy where possible)
# ---------------------------------------------------------------------------

C_LIGHT = 2.99792458e8            # m/s
HBAR = 1.054571817e-34            # J s
K_B = 1.380649e-23                # J/K
EV = 1.602176634e-19              # J
M_E_EV = 0.51099895000e6          # electron mass (eV)
M_P_EV = 938.27208816e6           # proton mass (eV)
M_N_EV = 939.56542052e6           # neutron mass (eV)
M_U_KG = 1.66053906660e-27        # atomic mass unit (kg)
M_P_KG = 1.67262192369e-27
M_E_KG = 9.1093837015e-31
SIGMA_T = 6.6524587321e-29        # Thomson cross-section (m^2)
R_E = 2.8179403262e-15            # classical electron radius (m)
ALPHA = 7.2973525693e-3
G_NEWTON = 6.67430e-11
MPC = 3.0856775814913673e22       # m
GPC = 1e3 * MPC
YEAR = 365.25 * 86400.0
ZETA3 = 1.2020569031595942
#: Reduced Planck mass and Planck mass in GeV.
M_PLANCK_GEV = 1.220890e19
#: Atomic masses (u) of 1H and 4He, for the mass per baryon (AME2020).
M_H1_U = 1.00782503223
M_HE4_U = 4.00260325413

# ---------------------------------------------------------------------------
# Cosmology we adopt (Planck 2018, TT,TE,EE+lowE+lensing; A&A 641, A6, Table 2)
# ---------------------------------------------------------------------------

#: CMB temperature today, Fixsen (2009, ApJ 707, 916).
T_CMB = 2.72548
T_CMB_ERR = 0.00057
H0 = 67.36                          # km/s/Mpc
OMEGA_M = 0.3153
OMBH2_PLANCK = 0.02237
OMBH2_PLANCK_ERR = 0.00015
N_EFF = 3.044
#: Primordial helium mass fraction for the eta <-> Omega_b h^2 conversion
#: (Planck's BBN-consistent value; the conversion moves by 1e-4 per 0.01).
YP = 0.2454


def h0_si(h0_kms_mpc: float = H0) -> float:
    return h0_kms_mpc * 1e3 / MPC


def omega_gamma_h2(t_cmb: float = T_CMB) -> float:
    """Photon density parameter times h^2."""
    a_rad = 4 * 5.670374419e-8 / C_LIGHT           # radiation constant (J m^-3 K^-4)
    rho_gamma = a_rad * t_cmb**4 / C_LIGHT**2       # kg/m^3
    rho_crit_h2 = 3 * (100e3 / MPC) ** 2 / (8 * math.pi * G_NEWTON)
    return rho_gamma / rho_crit_h2


def omega_r_h2(t_cmb: float = T_CMB, n_eff: float = N_EFF) -> float:
    """Radiation (photons + massless neutrinos) density times h^2."""
    return omega_gamma_h2(t_cmb) * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * n_eff)


# ---------------------------------------------------------------------------
# Figure styling
# ---------------------------------------------------------------------------

FIG_DPI = 200
