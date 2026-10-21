"""Project configuration: paths, unit conversions, the analysis choices, the
Galactic constants from sibling projects, and the published values we grade
ourselves against.

H4 measures the acceleration of the Solar System as it falls around the Milky
Way. Light from every source reaches us tilted by our own velocity
(aberration). As the Galaxy bends our path, that velocity changes, so the tilt
changes too, and every distant quasar appears to drift by ~5 microarcseconds
per year towards the direction of our acceleration -- the Galactic Centre.
Gaia measured the proper motions of 1.6 million quasars; the drift is a few
hundred times smaller than a single quasar's error bar, so only the average over
the whole sky reveals it.

Published answers live in the ``REFERENCE`` block at the bottom and are used
only to grade the results. Nothing in the analysis reads them.
"""

from __future__ import annotations

import math
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CHUNK_DIR = RAW_DIR / "chunks"
LIT_DIR = DATA_DIR / "literature"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, CHUNK_DIR, LIT_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Units
# --------------------------------------------------------------------------
#: Speed of light (m/s), the Julian year (s) and one microarcsecond (rad).
C_M_S = 299_792_458.0
JULIAN_YEAR_S = 365.25 * 86_400.0
UAS_RAD = math.pi / 180.0 / 3600.0 / 1e6
#: An acceleration a makes every quasar drift by g = a / c. One microarcsecond
#: per year of drift is therefore an acceleration of 4.606e-11 m/s^2
#: (Klioner et al. 2021 quote 4.60566e-11; test_physics.py checks it).
UAS_YR_TO_M_S2 = UAS_RAD / JULIAN_YEAR_S * C_M_S
#: ... or 1.453 km/s per million years.
UAS_YR_TO_KMS_MYR = UAS_YR_TO_M_S2 * 1e6 * JULIAN_YEAR_S / 1e3
#: A kiloparsec in metres (IAU 2012 au), and Newton's G in kpc (km/s)^2 / Msun.
KPC_M = 3.085677581491367e19
G_KPC_KMS2_MSUN = 4.300917270e-6
#: (km/s)^2 per kpc -> m/s^2 (an acceleration V^2/R written in Galactic units).
KMS2_PER_KPC_TO_M_S2 = 1e6 / KPC_M
#: 1 mas/yr at 1 kpc is 4.74047 km/s.
K_KMS = 4.740470463533348

# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------
#: ICRS -> Galactic rotation, the Hipparcos definition (ESA 1997, Vol. 1,
#: Sect. 1.5.3) that Gaia uses for its l, b columns: r_gal = A_G r_icrs.
#: test_physics.py checks it against the l, b that Gaia publishes.
A_G = (
    (-0.0548755604162154, -0.8734370902348850, -0.4838350155487132),
    (+0.4941094278755837, -0.4448296299600112, +0.7469822444972189),
    (-0.8676661490190047, -0.1980763734312015, +0.4559837761750669),
)
#: Sgr A* (Reid & Brunthaler 2004), ICRS degrees: the black hole at the
#: bottom of the Galaxy's potential, 3.3' and 2.7' from (l, b) = (0, 0).
SGRA_RA_DEG, SGRA_DEC_DEG = 266.4168371, -29.0078106

# --------------------------------------------------------------------------
# Data: Gaia DR3's quasar-like sources (the Gaia-CRF3)
# --------------------------------------------------------------------------
GAIA_SYNC_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
#: Sky chunks for the download: HEALPix level 1 (48 pixels, ~34,000 quasars
#: each); a chunk that times out is split into its four children, down to L4.
HEALPIX_BASE_LEVEL = 1
HEALPIX_MAX_LEVEL = 4
#: gaia_source columns fetched for every Gaia-CRF3 source. Astrometry and its
#: errors and correlations; the solution type and the colour the astrometry
#: used; quality indicators; photometry; Gaia's own classifier outputs.
GAIA_COLUMNS = [
    "source_id", "ra", "dec", "l", "b", "ecl_lat",
    "parallax", "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error",
    "pmra_pmdec_corr", "parallax_pmra_corr", "parallax_pmdec_corr",
    "astrometric_params_solved", "nu_eff_used_in_astrometry", "pseudocolour",
    "pseudocolour_error", "ruwe", "astrometric_excess_noise",
    "astrometric_excess_noise_sig", "astrometric_n_good_obs_al",
    "visibility_periods_used", "ipd_frac_multi_peak", "ipd_gof_harmonic_amplitude",
    "duplicated_source", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
    "bp_rp", "phot_bp_rp_excess_factor", "in_qso_candidates", "in_galaxy_candidates",
    "classprob_dsc_combmod_quasar", "classprob_dsc_combmod_galaxy",
    "classprob_dsc_combmod_star",
]
#: Columns from gaiadr3.agn_cross_id (which external catalogue put the source
#: in the Gaia-CRF3) and gaiadr3.qso_candidates (Gaia's own spectro-photometric
#: redshift, its quality flags, and whether a host galaxy was seen).
AGN_COLUMNS = ["catalogue_name"]
QSOC_COLUMNS = ["redshift_qsoc", "redshift_qsoc_lower", "redshift_qsoc_upper",
                "flags_qsoc", "host_galaxy_flag", "vari_best_class_name"]
#: Server-side counts the download must reproduce exactly (Klioner et al.
#: 2021, Sect. 4.1; checked again by fetch_data.py with COUNT(*) queries).
N_CRF3 = 1_614_173
N_CRF3_5P = 1_215_942
N_CRF3_6P = 398_231

#: SDSS DR16Q spectroscopic redshifts (Lyke et al. 2020, VizieR VII/289),
#: for the test that the drift does not depend on distance. Taken from project
#: I3's cache when it exists, otherwise fetched from VizieR in RA strips.
I3_DR16Q_PARQUET = ATLAS_DIR / "34-quasars" / "data" / "derived" / "dr16q_columns.parquet"
VIZIER_ASU_TSV = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
DR16Q_VIZIER = "VII/289/dr16q"
N_DR16Q = 750_414
#: Gaia positions (epoch 2016) and SDSS positions match to well under 1".
DR16Q_MATCH_ARCSEC = 1.0

# --------------------------------------------------------------------------
# The fit (Klioner et al. 2021's choices are the baseline; ours are variants)
# --------------------------------------------------------------------------
#: Highest degree of vector spherical harmonics fitted alongside the glide.
#: Klioner et al. found the glide stable from l_max = 3 and adopted it.
LMAX_BASELINE = 3
LMAX_SCAN = list(range(1, 11))
#: Outlier clipping: a source is rejected when its weighted post-fit residual
#: X exceeds KAPPA times the median X of all sources (Lindegren 2018).
KAPPA_BASELINE = 3.0
KAPPA_SCAN = [2.5, 3.0, 4.0, 5.0, 8.0, math.inf]
#: Bootstrap replicates for the uncertainties (each refits with clipping).
N_BOOTSTRAP = 400
#: Random seed for every resampling in the project.
SEED = 20261019
#: Magnitude bins of Klioner et al.'s Fig. 8 and their published counts.
G_BINS = [(-math.inf, 18.0), (18.0, 19.0), (19.0, 20.0), (20.0, math.inf)]
G_BIN_COUNTS_K21 = [29_200, 146_614, 490_161, 549_967]

# --------------------------------------------------------------------------
# The Galaxy: what the acceleration should be
# --------------------------------------------------------------------------
#: Project H1 (17-galactic-centre): the Galactic Centre's distance from the
#: orbit of S2, both teams' data, GR included -- 7.93 +- 0.04 (stat) +- 0.03
#: (sys) kpc. Read from H1's results.json when present (run.py checks).
R0_H1_KPC = 7.931
R0_H1_ERR_KPC = math.hypot(0.043, 0.03)
#: Project H2 (18-milky-way-rotation) did not measure the circular speed at
#: the Sun itself; it adopted Reid et al. (2019): R0 = 8.15 +- 0.15 kpc,
#: Theta0 = 236 +- 7 km/s, and the Sun's angular speed (Theta0 + V_sun)/R0 =
#: 30.32 +- 0.27 km/s/kpc.
R0_H2_KPC, R0_H2_ERR_KPC = 8.15, 0.15
THETA0_H2_KMS, THETA0_H2_ERR_KMS = 236.0, 7.0
OMEGA_SUN_H2, OMEGA_SUN_H2_ERR = 30.32, 0.27
#: Sgr A*'s apparent motion along the Galactic plane, the reflex of the Sun's
#: orbit (Reid & Brunthaler 2020): -6.411 +- 0.008 mas/yr. Times K this is
#: the Sun's angular speed around the Galactic Centre, 30.39 km/s/kpc.
MU_SGRA_L_MAS_YR, MU_SGRA_L_ERR = -6.411, 0.008
#: How much faster than a circular orbit the Sun moves along the rotation
#: (Schoenrich, Binney & Dehnen 2010): 12.24 +- 0.47 (stat) +- 2 (sys) km/s.
V_SUN_PEC_KMS, V_SUN_PEC_ERR_KMS = 12.24, 2.05
#: The Sun's height above the Galactic plane (Bland-Hawthorn & Gerhard 2016,
#: as used by Klioner et al.; Bennett & Bovy 2019 give 20.8 +- 0.3 pc).
Z_SUN_KPC, Z_SUN_ERR_KPC = 0.025, 0.005
#: The Large Magellanic Cloud: 49.59 kpc (Pietrzynski et al. 2019) towards
#: (l, b) = (280.47, -32.89) deg, total mass 1-2.5e11 Msun (Erkal et al. 2019).
LMC_D_KPC = 49.59
LMC_L_DEG, LMC_B_DEG = 280.47, -32.89
LMC_MASS_RANGE_MSUN = (1.0e11, 2.5e11)

# --------------------------------------------------------------------------
# REFERENCE: published values used only for grading
# --------------------------------------------------------------------------
REFERENCE = {
    "klioner2021": {
        "ref": "Gaia Collaboration, Klioner et al. 2021, A&A 649, A9 (arXiv:2012.02036)",
        # Table 2: bootstrap 1-sigma uncertainties, l_max = 3, kappa = 3.
        "g_icrs": (-0.07, -4.30, -2.64), "g_icrs_err": (0.41, 0.35, 0.36),
        "rho_icrs": {"xy": +0.001, "xz": -0.094, "yz": -0.025},
        "ra_dec": (269.1, -31.6), "ra_dec_err": (5.4, 4.1),
        "g_gal": (+5.04, -0.10, -0.29), "g_gal_err": (0.35, 0.36, 0.41),
        "rho_gal": {"XY": +0.036, "XZ": -0.014, "YZ": -0.079},
        "l_b": (358.9, -3.3), "l_b_err": (4.1, 4.6),
        "amplitude_uas_yr": (5.05, 0.35),
        "amplitude_kms_myr": (7.33, 0.51),
        "amplitude_1e10_m_s2": (2.32, 0.16),
        # Appendix A: |g| is biased high by about this much.
        "amplitude_bias_uas_yr": 0.024,
        # Sect. 6: unit weight error of all fits; outliers rejected at kappa=3.
        "uwe": 1.048, "n_outliers_approx": 3800,
        # Sect. 4.1 / Fig. 6: widths of the normalised parallax and proper
        # motion distributions (best-fit Gaussians).
        "normalised_widths": {"parallax": 1.052, "pmra": 1.055, "pmdec": 1.063},
        # Table 1: medians for the five-parameter sources.
        "table1_5p": {"G": 19.92, "bp_rp": 0.64, "nu_eff": 1.589, "ruwe": 1.013,
                      "pmra_error_uas": 457, "pmdec_error_uas": 423},
        "table1_6p": {"G": 20.46, "bp_rp": 0.92, "ruwe": 1.044,
                      "pmra_error_uas": 892, "pmdec_error_uas": 832},
        # Sect. 3: the expected acceleration, km/s/Myr, galactic (X, Y, Z):
        # centripetal from R0 = 8.178 kpc, Sgr A*'s proper motion and
        # Schoenrich's V_sun; bar and spirals; pull towards the plane.
        "expected_kms_myr": (+6.98, -0.06, -0.13), "expected_err": (0.12, 0.05, 0.03),
        "expected_centripetal_uas_yr": (4.80, 0.08),
        "expected_with_lmc_range_kms_myr": ((+6.89, -0.20, -0.20), (+7.17, -0.48, -0.40)),
    },
    # VLBI: the ICRF3 data set (Charlot et al. 2020, A&A 644, A159).
    "charlot2020": {"ref": "Charlot et al. 2020, A&A 644, A159 (ICRF3)",
                    "amplitude_uas_yr": (5.83, 0.23), "ra_dec": (270.2, -20.2),
                    "ra_dec_err": (2.3, 3.6)},
    "titov_krasna2018": {"ref": "Titov & Krasna 2018, A&A 610, A36",
                         "amplitude_uas_yr": (5.2, 0.2)},
    # Re-analyses of the same Gaia EDR3/DR3 astrometry, 2025-26 (values checked
    # in the papers' arXiv sources on 19 Oct 2026).
    "tsigkas_kouvelis2026": {"ref": "Tsigkas Kouvelis et al. 2026, JCAP 06, 034 (arXiv:2508.02810)",
                             "amplitude_uas_yr": (4.86, 0.34), "ra_dec": (266.9, -31.6),
                             "note": "Bayesian VSH to l = 4; drift ~2x larger at z > 1.25 than below (2-3 sigma)"},
    "nguyen_murray2026": {"ref": "Nguyen & Murray 2026, A&A in press (arXiv:2605.30367)",
                          "amplitude_uas_yr": (4.94, 0.32), "ra_dec": (267.9, -34.4),
                          "note": "simulation-based inference, l_max = 3, EDR3 five-parameter sources"},
    "darling2025": {"ref": "Darling 2025, ApJL 982, L46 (arXiv:2412.08605)",
                    "amplitude_uas_yr": (5.20, 0.31), "ra_dec": (273.8, -15.8),
                    "note": "Quaia five-parameter quasars, dipole (l = 1) only"},
    "makarov2025": {"ref": "Makarov 2025, Nat. Astron. 9, 1396 (arXiv:2507.01798)",
                    "gc_component_uas_yr": (5.39, 0.38), "note": "Galactocentric component, quasars at z > 1"},
    "macmillan2019": {"ref": "MacMillan et al. 2019, A&A 630, A93 (IVS working group)",
                      "recommended_uas_yr": 5.8},
    "gravity2019_r0": {"ref": "GRAVITY Collaboration 2019, A&A 625, L10",
                       "r0_kpc": (8.178, math.hypot(0.013, 0.022))},
    "gravity2022_r0": {"ref": "GRAVITY Collaboration 2022, A&A 657, L12",
                       "r0_kpc": (8.277, math.hypot(0.009, 0.033))},
    "eilers2019": {"ref": "Eilers et al. 2019, ApJ 871, 120",
                   "vc_kms": (229.0, 0.2), "vc_sys_frac": 0.02, "r0_kpc": 8.122},
}
