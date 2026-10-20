"""Project configuration: paths, physical constants, inputs from earlier projects,
and the published numbers we grade ourselves against.

N6 downloads almost nothing. Its data are the results of earlier projects in
this atlas, read from their own folders:

* H1 (17-galactic-centre): the S2 measurements and the posterior chain of the
  joint Keck + VLT relativistic orbit fit;
* I5 (36-m87-shadow): the measured ring diameters of M87* and Sgr A*;
* I4 (35-x-ray-binaries): the masses of the stellar black holes;
* D7 (37-expansion-tests): the 46 published measurements of the CMB temperature
  at redshifts up to 6.3;
* C1, I1-I3 (02-gravitational-waves, 32-, 33-, 34-): more black-hole masses.

The one new download is the LIGO-Virgo-KAGRA posterior samples for GW250114,
the loudest black-hole merger so far (27 MB), standing in for project I7,
which has not been built yet. ``fetch_data.py`` records the SHA-256 of every
input, so a later change in an earlier project's results is noticed.

If an earlier project's results are missing, the analysis falls back to that
project's published headline numbers (listed below) and says so.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"
VIEWER_PATH = FIGURE_DIR / "kerr_orbit_explorer.html"

for _d in (DATA_DIR, RAW_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical constants: CODATA 2018 and the IAU 2015 nominal solar values
# ---------------------------------------------------------------------------

C = 299_792_458.0                 # speed of light, m/s (exact)
G = 6.674_30e-11                  # Newton's constant, m^3 kg^-1 s^-2
HBAR = 1.054_571_817e-34          # reduced Planck constant, J s (exact)
K_B = 1.380_649e-23               # Boltzmann constant, J/K (exact)
SIGMA_SB = 5.670_374_419e-8       # Stefan-Boltzmann constant, W m^-2 K^-4
GM_SUN = 1.327_124_400_18e20      # heliocentric gravitational constant, m^3/s^2 (IAU 2015)
M_SUN = GM_SUN / G                # 1.98841e30 kg
AU = 149_597_870_700.0            # astronomical unit, m (IAU 2012, exact)
PC = AU * 648_000.0 / 3.141592653589793
YEAR = 365.25 * 86_400.0          # Julian year, s
GYR = 1e9 * YEAR
L_PLANCK2 = G * HBAR / C**3       # Planck length squared, m^2
RAD_TO_UAS = 180.0 / 3.141592653589793 * 3600.0 * 1e6

#: Today's CMB temperature (Fixsen 2009, ApJ 707, 916).
T_CMB = 2.72548

#: Planck 2018 base-LCDM (TT,TE,EE+lowE+lensing; Planck Collaboration 2020,
#: A&A 641, A6, Table 2 column 5): used only for illustrations of the
#: expanding universe (photon redshift, peculiar velocities, the far future).
PLANCK18 = {"H0": 67.36, "Omega_m": 0.3153, "Omega_r": 9.1e-5}

# ---------------------------------------------------------------------------
# Inputs from earlier projects
# ---------------------------------------------------------------------------

H1_DIR = ATLAS_DIR / "17-galactic-centre"
H1_RESULTS = H1_DIR / "data" / "results.json"
H1_CHAIN = H1_DIR / "data" / "chains" / "joint_gr.npz"
I5_RESULTS = ATLAS_DIR / "36-m87-shadow" / "data" / "results.json"
I4_RESULTS = ATLAS_DIR / "35-x-ray-binaries" / "data" / "results.json"
D7_TCMB_TABLE = ATLAS_DIR / "37-expansion-tests" / "literature" / "tcmb_measurements.csv"
D7_RESULTS = ATLAS_DIR / "37-expansion-tests" / "data" / "results.json"
C1_RESULTS = ATLAS_DIR / "02-gravitational-waves" / "data" / "results.json"
I1_RESULTS = ATLAS_DIR / "32-m-sigma-relation" / "data" / "results.json"
I2_RESULTS = ATLAS_DIR / "33-reverberation-mapping" / "data" / "results.json"
I3_RESULTS = ATLAS_DIR / "34-quasars" / "data" / "derived" / "results.json"
#: Wu & Shen (2022) virial masses of the SDSS DR16 quasars, as cached by I3.
I3_WUSHEN = ATLAS_DIR / "34-quasars" / "data" / "derived" / "wushen22_columns.parquet"
#: Kormendy & Ho (2013) black-hole masses, as parsed and cached by I1.
I1_KH13 = ATLAS_DIR / "32-m-sigma-relation" / "data" / "cache" / "kh13.csv"

INPUT_FILES = {
    "H1 results": H1_RESULTS, "H1 joint GR chain": H1_CHAIN,
    "H1 Keck astrometry": H1_DIR / "data" / "raw" / "do2019" / "S0-2_astro.csv",
    "H1 Keck velocities": H1_DIR / "data" / "raw" / "do2019" / "S0-2_rv.csv",
    "H1 VLT table": H1_DIR / "data" / "raw" / "gillessen2017_table5.dat",
    "I5 results": I5_RESULTS, "I4 results": I4_RESULTS,
    "D7 CMB temperature table": D7_TCMB_TABLE, "D7 results": D7_RESULTS,
    "C1 results": C1_RESULTS, "I1 results": I1_RESULTS, "I2 results": I2_RESULTS,
    "I3 results": I3_RESULTS, "I3 Wu & Shen masses": I3_WUSHEN, "I1 Kormendy & Ho table": I1_KH13,
}

# ---------------------------------------------------------------------------
# The one download: GW250114 posterior samples (LIGO-Virgo-KAGRA)
# ---------------------------------------------------------------------------

GW250114 = {
    "name": "GW250114_082203",
    "paper": ("LIGO Scientific, Virgo and KAGRA Collaborations (Abac, A. G. et al.) 2025, "
              "Phys. Rev. Lett. 135, 111403, 'GW250114: Testing Hawking's Area Law and the "
              "Kerr Nature of Black Holes' (arXiv:2509.08054; doi:10.1103/kw5g-d732)"),
    "zenodo_record": "https://doi.org/10.5281/zenodo.16877102",
    # The discovery paper's preferred analysis (GWOSC flags it 'preferred').
    "samples_file": "posterior_samples_NRSur7dq4.h5",
    "samples_url": ("https://zenodo.org/api/records/16877102/files/"
                    "posterior_samples_NRSur7dq4.h5/content"),
    "gwosc_parameters_url": ("https://gwosc.org/api/v2/event-versions/"
                             "GW250114_082203-v1/parameters"),
}

USER_AGENT = "astro-hobby-project/1.0 (63-killing-vectors; python-urllib)"

# ---------------------------------------------------------------------------
# Fallbacks: earlier projects' headline numbers, used only if a file is missing
# ---------------------------------------------------------------------------

FALLBACK = {
    # H1: joint Keck+VLT GR fit (17-galactic-centre/README.md)
    "sgra_mass_1e6": (3.957, 0.042), "sgra_r0_kpc": (7.931, 0.043),
    # I5 (36-m87-shadow/README.md): ring diameters, micro-arcseconds
    "m87_ring_uas": (43.8, 0.4), "sgra_ring_uas": (53.9, 1.3),
}

# ---------------------------------------------------------------------------
# Textbook results we must reproduce (G = c = 1, M = 1 unless stated)
# ---------------------------------------------------------------------------

#: Number of independent Killing vectors. Flat, Schwarzschild, Kerr and FLRW
#: are in every GR textbook (e.g. Carroll 2004, Spacetime and Geometry,
#: sections 3.8, 5.1, 6.6 and 8.1). The others are validation cases with
#: well-known answers: Goedel's rotating universe has a 5-dimensional isometry
#: group (Goedel 1949, Rev. Mod. Phys. 21, 447), the Einstein static universe
#: R x S^3 a 7-dimensional one, and a generic plane gravitational wave 5
#: (Bondi, Pirani & Robinson 1959, Proc. R. Soc. A 251, 519). De Sitter space
#: is maximally symmetric (10) even in its expanding FLRW form.
KILLING_COUNTS = {
    "Minkowski": 10, "Minkowski (spherical coordinates)": 10, "Schwarzschild": 4,
    "Kerr": 2, "FLRW, generic a(t)": 6, "FLRW, matter era": 6, "FLRW, radiation era": 6,
    "FLRW, Planck 2018 LCDM": 6, "de Sitter (FLRW form)": 10, "Goedel universe": 5,
    "Einstein static universe": 7, "plane gravitational wave": 5,
    # static and axially symmetric only: d_t and d_phi (the metric depends on neither t nor phi)
    "Zipoy-Voorhees (delta = 2)": 2,
}

PUBLISHED = {
    "r_isco_schwarzschild": 6.0,
    "r_photon_schwarzschild": 3.0,
    "b_crit_schwarzschild": 27 ** 0.5,
    "efficiency_schwarzschild": 1 - (8 / 9) ** 0.5,        # 5.72%
    "efficiency_extremal_prograde": 1 - 1 / 3 ** 0.5,       # 42.3%
    "efficiency_extremal_retrograde": 1 - (25 / 27) ** 0.5,  # 3.78%
    "r_isco_extremal_retrograde": 9.0,
    "max_rotational_energy_fraction": 1 - 2 ** -0.5,        # 29.3%
    "hawking_T_1msun_K": 6.17e-8,
}
