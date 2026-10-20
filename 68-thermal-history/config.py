"""Shared settings: units, constants, data locations and the literature values we check against.

Units. Inside the code every energy, mass and temperature is in **GeV** (natural units,
k_B = hbar = c = 1), so a temperature of 1 MeV is written 1e-3. Times come out in GeV^-1
and are converted to seconds with hbar. 1 MeV of temperature is 1.160 x 10^10 K.
"""

from __future__ import annotations

import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
RAW_DIR = DATA_DIR / "raw"              # tables downloaded as they are
LIT_DIR = DATA_DIR / "literature"       # arXiv sources (tables read straight from the TeX)
DERIVED_DIR = DATA_DIR / "derived"      # our own computed curves (npz)
FIG_DIR = HERE / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

ROOT = HERE.parent  # the atlas; sibling projects' results are reused when present

# --- Constants (CODATA 2022 / PDG 2024) ----------------------------------------------

HBAR_GEV_S = 6.582119569e-25          # hbar in GeV s
K_B_GEV_PER_K = 8.617333262e-14       # Boltzmann constant in GeV/K (exact)
C_M_S = 299_792_458.0                 # m/s (exact)
HBARC_GEV_M = 1.973269804e-16         # hbar c in GeV m
G_SI = 6.67430e-11                    # m^3 kg^-1 s^-2 (CODATA, +- 0.00015)
GEV_KG = 1.78266192e-27               # 1 GeV/c^2 in kg
HBAR_J_S = 1.054571817e-34           # hbar in J s (exact)
#: Planck mass sqrt(hbar c / G) in GeV = 1.220890 x 10^19 GeV (PDG: 1.220890(14) x 10^19).
M_PL = math.sqrt(HBAR_J_S * C_M_S / G_SI) / GEV_KG
M_SUN_KG = 1.988416e30                # IAU 2015 nominal GM_sun / G (CODATA 2018 G)
YEAR_S = 365.25 * 86_400.0            # Julian year
MPC_M = 3.0856775814913673e22
ALPHA_EM = 1 / 137.035999177          # fine-structure constant (CODATA 2022)
G_FERMI = 1.1663788e-5                # Fermi constant, GeV^-2 (PDG)
M_E = 0.51099895069e-3                # electron mass, GeV (CODATA 2022)

ZETA3 = 1.2020569031595942

# --- Today's universe ---------------------------------------------------------------

#: CMB temperature today: Fixsen (2009, ApJ 707, 916), from FIRAS recalibrated with WMAP.
T0_K = 2.72548
T0_K_ERR = 0.00057
T0_GEV = T0_K * K_B_GEV_PER_K

#: Late-time cosmology for the calendar after nucleosynthesis: Planck 2018
#: TT,TE,EE+lowE+lensing+BAO (Planck Collaboration 2020, A&A 641, A6, Table 2, last column),
#: the same universe as astropy's Planck18 and the atlas's O1/O2 projects.
PLANCK18 = {
    "H0": 67.66,            # km/s/Mpc
    "ombh2": 0.02242,
    "omch2": 0.11933,
    "m_nu_eV": 0.06,        # one massive neutrino (the rest massless), Planck's baseline
    "source": "Planck 2018 VI, Table 2, TT,TE,EE+lowE+lensing+BAO (= astropy Planck18)",
}

# --- Temperature range of the thermal history ----------------------------------------

T_MAX = 1.0e4      # 10 TeV: the top of the staircase (GeV)
T_MIN = 1.0e-5     # 10 keV: e+e- annihilation is over (GeV)
#: Grid points per decade for all curves.
N_PER_DECADE = 240

# --- QCD -----------------------------------------------------------------------------

#: Hadron-resonance gas: keep hadrons below this mass (S&S used < 2.5 GeV; heavier states add
#: < 1e-4 to g* below 200 MeV).
HRG_MASS_MAX = 3.0
#: f0(500) ("sigma") and K0*(700) ("kappa") are broad, and their attraction is cancelled by
#: repulsive channels (Broniowski, Giacosa & Begun 2015, PRC 92, 034905); hadron-gas models
#: usually leave them out. We do too, and report the effect of keeping them.
HRG_EXCLUDE = (9000221, 9000311, -9000311, 9000321, -9000321)
#: K_S and K_L are mixtures of K0 and anti-K0, which the PDG table also lists: drop them.
HRG_DUPLICATES = (130, 310)

#: Hadron gas below this temperature, lattice QCD above (S&S and Borsanyi et al. 2016: 120 MeV).
T_HRG_TO_LATTICE = 0.120
#: Lattice QCD below this temperature, perturbative QCD above (Borsanyi et al. 2016 matched at
#: 500 MeV - 1 GeV; S&S tried both).
T_LATTICE_TO_PQCD = 1.0
#: Strong coupling at the Z mass (PDG 2024 world average 0.1180 +- 0.0009).
ALPHA_S_MZ = 0.1180
M_Z = 91.1880
#: The unknown g^6 constant of the QCD pressure, fixed by Borsanyi et al. 2016 against their
#: lattice data (also S&S's choice).
QC_PQCD = -3000.0

#: The textbook, ideal-gas picture switches abruptly from free quarks and gluons to hadrons
#: at the lattice crossover temperature.
T_SWITCH_IDEAL = 0.1565

# --- Neutrino decoupling --------------------------------------------------------------

#: Effective low-energy neutrino-electron couplings including radiative corrections, as used by
#: Escudero (2020, JCAP 05, 048, eq. 3.6): g_eL, g_eR for nu_e; g_muL, g_muR for nu_mu,tau.
G_EL, G_ER, G_MUL, G_MUR = 0.727, 0.233, -0.273, 0.233
#: Fermi-Dirac corrections to the Maxwell-Boltzmann energy-transfer rates (Escudero 2020, eq. 3.7).
F_A_FD, F_S_FD = 0.884, 0.829
#: Start of the decoupling integration (all species in equilibrium), and end.
T_NU_START = 10e-3
T_NU_END = 5e-6

# --- Literature values ---------------------------------------------------------------

LIT = {
    "g_high_ideal": {"value": 106.75, "source": "Standard Model ideal gas, all species relativistic (spec)"},
    "g_few_MeV": {"value": 10.75, "source": "photons + e+- + 3 neutrinos (spec)"},
    "g_rho_today_inst": {"value": 3.36, "source": "instantaneous decoupling, N_eff = 3 (spec)"},
    "g_s_today_inst": {"value": 3.91, "source": "instantaneous decoupling (spec)"},
    "T_pc_hotqcd2019": {"value": 0.1565, "err": 0.0015,
                        "source": "HotQCD (Bazavov et al. 2019, PLB 795, 15; arXiv:1812.08235): chiral crossover"},
    "T_c_borsanyi2020": {"value": 0.1580, "err": 0.0006, "width": 0.015, "width_err": 0.001,
                         "source": "Borsanyi et al. 2020 (PRL 125, 052001; arXiv:2002.02821)"},
    "T_ew_donofrio2016": {"value": 159.5, "err": 1.5,
                          "source": "D'Onofrio & Rummukainen 2016 (PRD 93, 025003; arXiv:1508.07161)"},
    "N_eff_bennett2021": {"value": 3.0440, "err": 0.0002,
                          "source": "Bennett et al. 2021 (JCAP 04, 073; arXiv:2012.02726)"},
    "N_eff_escudero2020": {"value": 3.045, "common_nu": 3.045, "separate_nu": 3.044,
                           "fd_nlo_common": 3.049, "fd_nlo_separate": 3.048,
                           "ratio_fd_nlo_common": 1.3954, "ratio_fd_nlo_separate_e": 1.3941,
                           "ratio_fd_nlo_separate_mu": 1.3962,
                           "g_s0": 3.931, "g_rho0": 3.383,
                           "source": "Escudero 2020 (JCAP 05, 048; arXiv:2001.04466), Table 1"},
    "T_nu_today_K": {"value": 1.945, "source": "(4/11)^(1/3) x 2.7255 K (spec)"},
    "n_nu_today_cm3": {"value": 336, "source": "112 per flavour, nu + anti-nu (P3 spec)"},
}

#: Annotated marks for the calendar. The energies of the quark-gluon plasma made in the lab are
#: "effective temperatures" from thermal photons (blue-shifted by the plasma's flow).
LAB_MARKS = {
    # filled in from a primary-source check; see README "Data credits"
}
