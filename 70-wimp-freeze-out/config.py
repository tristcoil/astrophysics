"""Shared settings: units, constants, data locations and the literature values we check against.

Units. Inside the code every energy, mass and temperature is in **GeV** (natural units,
k_B = hbar = c = 1), so a temperature of 1 MeV is written 1e-3. Cross sections times
velocity are carried in GeV^-2 and converted to cm^3/s only for input and output.
"""

from __future__ import annotations

import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
RAW_DIR = DATA_DIR / "raw"              # files downloaded as they are
LIT_DIR = DATA_DIR / "literature"       # members extracted from arXiv sources
DERIVED_DIR = DATA_DIR / "derived"      # our own computed curves
FIG_DIR = HERE / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

ROOT = HERE.parent  # the atlas; sibling projects' caches are reused when present

# --- Constants (CODATA 2022 / PDG 2024) ---------------------------------------------------

HBAR_GEV_S = 6.582119569e-25          # hbar in GeV s
HBARC_GEV_CM = 1.973269804e-14        # hbar c in GeV cm
C_CM_S = 2.99792458e10                # cm/s (exact)
K_B_GEV_PER_K = 8.617333262e-14       # Boltzmann constant in GeV/K (exact)
G_SI = 6.67430e-11                    # m^3 kg^-1 s^-2
HBAR_J_S = 1.054571817e-34            # J s (exact)
C_M_S = 299_792_458.0
GEV_KG = 1.78266192e-27               # 1 GeV/c^2 in kg
#: Planck mass sqrt(hbar c / G) = 1.220890e19 GeV (PDG 1.220890(14)e19).
M_PL = math.sqrt(HBAR_J_S * C_M_S / G_SI) / GEV_KG
#: Reduced Planck mass M_Pl / sqrt(8 pi) = 2.435e18 GeV (Saikawa & Shirai 2020 use 2.4353e18).
M_PL_RED = M_PL / math.sqrt(8 * math.pi)
#: <sigma v> conversion: 1 GeV^-2 = (hbar c)^2 c = 1.1673e-17 cm^3/s.
GEV2_TO_CM3S = HBARC_GEV_CM ** 2 * C_CM_S
#: Cross-section conversion: 1 GeV^-2 = (hbar c)^2 = 3.894e-28 cm^2.
GEV2_TO_CM2 = HBARC_GEV_CM ** 2
ZETA3 = 1.2020569031595942

# --- Today's universe -----------------------------------------------------------------

#: CMB temperature: Fixsen (2009, ApJ 707, 916). Saikawa & Shirai and Planck use 2.7255 K.
T0_K = 2.7255
T0_GEV = T0_K * K_B_GEV_PER_K
#: Critical density / h^2 = 1.05371255(23) x 10^-5 GeV cm^-3 (PDG 2024, Astrophysical constants).
RHO_CRIT_H2_GEV_CM3 = 1.05371255e-5

#: Planck 2018 cold dark matter density (Planck Collaboration 2020, A&A 641, A6, Table 2,
#: TT,TE,EE+lowE+lensing). Saikawa & Shirai 2020 used 0.120 +- 0.001.
OMEGA_C_H2 = 0.1200
OMEGA_C_H2_ERR = 0.0012
#: Steigman, Dasgupta & Beacom (2012) normalised to WMAP-7: Omega h^2 = 0.11.
OMEGA_H2_SDB12 = 0.11

# --- Freeze-out solver ----------------------------------------------------------------

#: Start the integration at x = m/T = X_START, deep in equilibrium.
X_START = 3.0
#: ... and stop at T = 1 eV (CMB era), as Saikawa & Shirai (2020) do.
T_END = 1.0e-9
#: Relative tolerance of the stiff solver (on ln Y).
RTOL = 1e-9
#: Masses for the thermal cross-section curve (GeV). The spec asks for 1 GeV - 10 TeV;
#: we go further on both sides to show where the curve flattens and where the
#: Saikawa & Shirai table ends.
MASS_MIN, MASS_MAX, N_MASS = 0.01, 1.0e5, 141   # 20 per decade
#: The particle types of Saikawa & Shirai (2020): internal degrees of freedom g.
DM_TYPES = {
    "real_scalar": {"g": 1, "label": "real scalar", "ss2020": "S_1"},
    "majorana": {"g": 2, "label": "Majorana fermion / complex scalar", "ss2020": "S_2"},
    "vector": {"g": 3, "label": "real vector", "ss2020": "S_3"},
    "dirac": {"g": 4, "label": "Dirac fermion", "ss2020": "S_4"},
}

# --- Neutron freeze-out (the same physics at ~1 MeV) -----------------------------------

M_E = 0.51099895069e-3                # electron mass, GeV (CODATA 2022)
M_N = 0.93956542194                   # neutron mass, GeV (CODATA 2022)
M_P = 0.93827208943                   # proton mass, GeV (CODATA 2022)
Q_NP = M_N - M_P                      # 1.29333 MeV
M_U = 0.93149410372                   # atomic mass unit, GeV (CODATA 2022)
#: Neutron lifetime: PDG 2024 average 878.4 +- 0.5 s (bottle experiments dominate).
TAU_N = 878.4
TAU_N_ERR = 0.5
#: The beam-method value (Yue et al. 2013, PRL 111, 222501): the "neutron lifetime puzzle".
TAU_N_BEAM = 887.7
#: Baryon density for the nucleosynthesis run: Planck 2018 TT,TE,EE+lowE+lensing.
OMEGA_B_H2 = 0.02237
#: Neutrino species and their effective number (Bennett et al. 2021: 3.0440).
N_EFF = 3.044

# --- Gamma rays from dwarf galaxies -------------------------------------------------------

#: McDaniel et al. (2024): TS profiles on 40 masses (1 GeV - 10 TeV) x 60 cross sections
#: (1e-28 - 1e-22 cm^3/s); the upper limit is where TS falls 2.71 below its maximum.
MCD_MASSES = (1.0, 1.0e4, 40)
MCD_SIGMAV = (1.0e-28, 1.0e-22, 60)
DELTA_TS_95 = 2.71
#: Their samples (notebook of the data release).
DWARF_SAMPLES = ("benchmark", "measured")
#: PPPC 4 DM ID tables start at 5 GeV.
GAMMA_MASS_MIN = 5.0

# --- Literature values we check against --------------------------------------------------

LIT = {
    "sdb12_sigmav_highmass": {"value": 2.2e-26, "source": "Steigman, Dasgupta & Beacom 2012 (PRD 86, "
                              "023506): 2.2e-26 cm^3/s above ~10 GeV for Omega h^2 = 0.11"},
    "sdb12_sigmav_peak": {"value": 5.2e-26, "mass": 0.3, "source": "SDB12 abstract: maximum near 0.3 GeV"},
    "canonical_sigmav": {"value": 3.0e-26, "source": "the textbook 'canonical' value"},
    "ss2020_sigmav_100gev": {"value": 2.11e-26, "source": "Saikawa & Shirai 2020 (JCAP 08, 011), "
                             "Majorana, s-wave, Omega h^2 = 0.120"},
    "ackermann15_crossing_gev": {"value": 100.0, "source": "Ackermann et al. 2015 (PRL 115, 231301): "
                                 "limits below the thermal cross section for masses <~ 100 GeV"},
    "mcdaniel24_excess": {"bb_mass": (150, 230), "tau_mass": (30, 50), "local_sigma": 2.0,
                          "global_sigma_bb": 0.5, "source": "McDaniel et al. 2024 (PRD 109, 063024)"},
    "lz2024_min": {"value": 2.2e-48, "mass": 40.0, "source": "LZ 2024 (PRL 135, 011802): 90% CL"},
    "yp_aver2021": {"value": 0.2453, "err": 0.0034, "source": "Aver et al. 2021 (JCAP 03, 027): "
                    "helium in metal-poor galaxies"},
    "br_inv_atlas2023": {"value": 0.107, "source": "ATLAS 2023 (PLB 842, 137963): BR(H -> invisible) < 0.107, "
                         "95% CL, Run 1 + 2"},
    "dh_cooke2018": {"value": 2.527e-5, "err": 0.030e-5, "source": "Cooke, Pettini & Steidel 2018 "
                     "(ApJ 855, 102)"},
}
