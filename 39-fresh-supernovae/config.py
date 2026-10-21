"""Project configuration: paths, data sources, sample cuts, and reference values.

Everything tunable lives here, together with the published numbers we compare
against. Published values are used only for comparison, never as fit inputs.
The one deliberate exception is D1: this project's whole point is to test D1's
standardisation (alpha, beta, M) and cosmology on supernovae D1 never saw.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
ZTF_DIR = DATA_DIR / "ztf_dr2"
ZTF_LC_DIR = ZTF_DIR / "lightcurves"
PANTHEON_DIR = DATA_DIR / "pantheon"
DUST_DIR = DATA_DIR / "sfd"
SNCOSMO_DIR = DATA_DIR / "sncosmo"
LIVE_DIR = DATA_DIR / "live"
CACHE_DIR = DATA_DIR / "cache"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR, ZTF_LC_DIR, PANTHEON_DIR, DUST_DIR, SNCOSMO_DIR,
           LIVE_DIR, CACHE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: Sibling project D1, whose Hubble diagram and cosmology we test.
D1_DIR = PROJECT_DIR.parent / "05-supernova-cosmology"

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: ZTF SN Ia DR2 (Rigault et al. 2025, A&A 694, A1): one zip with every light
#: curve, spectrum and table. We read only the members we need, by byte range.
ZTF_DR2_ZIP = "https://ztfcosmo.in2p3.fr/download/data"
ZTF_DR2_PREFIX = "ztfsniadr2_lite/"

#: Tables taken from the release (paths inside the zip, after the prefix).
ZTF_DR2_TABLES = [
    "README.md",
    # SALT2-T21 fits by the ZTF team (their default), redshifts, sub-types,
    # quality flags, Milky Way E(B-V).
    "tables/snia_data.csv",
    # Host galaxy stellar mass and colour: whole galaxy, and within 2 kpc.
    "tables/globalhost_data.csv",
    "tables/localhost_data.csv",
    # The team's own SALT3 fits, same phase range and bands as ours: the
    # answer key for our light-curve fitting.
    "tables/extra/snia_data_phase-10to40_gri_salt3.csv",
    # Their SALT2 fits with only g and r, to see what the i band adds.
    "tables/extra/snia_data_phase-10to40_gr_salt2-T21.csv",
]

#: Pantheon+ (Scolnic et al. 2022; Brout et al. 2022), the sample D1 fitted.
#: If D1's cached copy exists and its SHA-256 matches, it is copied instead.
PANTHEON_RAW = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
                "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")
PANTHEON_FILES = {
    "Pantheon+SH0ES.dat": PANTHEON_RAW + "Pantheon%2BSH0ES.dat",
    "Pantheon+SH0ES_STAT+SYS.cov": PANTHEON_RAW + "Pantheon%2BSH0ES_STAT%2BSYS.cov",
}

#: Schlegel, Finkbeiner & Davis (1998) dust maps, in the format the `sfdmap`
#: package reads (north and south Galactic polar projections, 4096^2 pixels).
SFD_RAW = "https://github.com/kbarbary/sfddata/raw/master/"
SFD_FILES = ["SFD_dust_4096_ngp.fits", "SFD_dust_4096_sgp.fits"]

#: Schlafly & Finkbeiner (2011) recalibration of SFD: E(B-V)_true = 0.86 x SFD.
#: Pantheon+ uses the same scaling with the Fitzpatrick (1999) law and R_V = 3.1.
SFD_SCALING = 0.86
MW_RV = 3.1

#: D1's light-curve model: the SALT2 surfaces retrained for Pantheon+ (Brout et
#: al. 2022, "B21"/"B22"), from the Pantheon+ release (3_SALT2, MODEL000).
SALT2_B21_DIR = PANTHEON_DIR / "SALT2_B21"
SALT2_B21_RAW = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
                 "Pantheon%2B_Data/3_SALT2/SALT2_B21trained_withsys/SALT2.MODEL000/")
SALT2_B21_FILES = ["SALT2.INFO", "salt2_color_correction.dat.gz",
                   "salt2_color_dispersion.dat.gz", "salt2_lc_dispersion_scaling.dat.gz",
                   "salt2_lc_relative_covariance_01.dat.gz", "salt2_lc_relative_variance_0.dat.gz",
                   "salt2_lc_relative_variance_1.dat.gz", "salt2_template_0.dat.gz",
                   "salt2_template_1.dat.gz"]

#: A sample of Pantheon+'s own Foundation light curves (PS1 griz, SNANA
#: format), refitted here to measure how sncosmo's amplitudes relate to
#: SNANA's (the MAG_OFFSET check). Fetched from the release's photometry folder.
FOUNDATION_LC_DIR = PANTHEON_DIR / "foundation_lcs"
FOUNDATION_LC_RAW = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
                     "Pantheon%2B_Data/1_DATA/photometry/Foundation_DJ17/")
N_FOUNDATION_CHECK = 40

#: SALT model and passbands, downloaded by sncosmo into data/sncosmo on first use.
SALT_SOURCE = "salt3"  # Kenworthy et al. 2021, the model the spec asks for
SALT_VERSION = "2.0"
ZTF_BANDS = {"ztfg": "ztfg", "ztfr": "ztfr", "ztfi": "ztfi"}

# --- The live element: supernovae classified in the last few weeks ---------

#: Transient Name Server public search (no account needed for this CSV view).
#: objtype 3 = "SN Ia". Only spectroscopically classified objects are returned.
TNS_SEARCH = "https://www.wis-tns.org/search"
TNS_OBJTYPE_SN_IA = 3

#: ALeRCE broker: object search, alert detections and the forced photometry
#: that ZTF now ships inside its alerts (30 days of history per alert).
ALERCE_API = "https://api.alerce.online/ztf/v1"
ALERCE_FORCED = "https://api.alerce.online/v2/lightcurve/forced-photometry/{oid}"

#: How far back the live snapshot looks, by discovery date. TNS's public search
#: filters on discovery date, not classification date, so "classified in the
#: last month" is taken as "discovered in the last 45 days and classified".
LIVE_DISCOVERY_WINDOW_DAYS = 45

# ---------------------------------------------------------------------------
# D1's fitted universe (05-supernova-cosmology, flat LambdaCDM with the SH0ES
# Cepheid calibration). Read from D1's results.json when it exists; these
# copies of D1's best fit are the fallback if D1 has not been run.
# ---------------------------------------------------------------------------

D1_FALLBACK = {"Om": 0.332, "H0": 73.531, "M": -19.24404}

# ---------------------------------------------------------------------------
# Sample cuts: the ZTF DR2 "basic cuts" (Rigault et al. 2025, Table 1),
# applied to our own SALT3 fits.
# ---------------------------------------------------------------------------

CUT_X1 = (-3.0, 3.0)
CUT_C = (-0.2, 0.8)
CUT_T0_ERR = 1.0
CUT_X1_ERR = 1.0
CUT_C_ERR = 0.1
CUT_FITPROB = 1e-7

#: Redshift windows used for statistics.
#: - volume-limited: below z = 0.06 ZTF found essentially every normal SN Ia
#:   it could classify, so there is no Malmquist bias (Amenouche et al. 2025);
#: - the "Hubble flow": above z = 0.023 peculiar velocities are < 7% of the
#:   recession velocity (the SH0ES window starts here, as in D1);
#: - DR2's own clean-scatter window, 0.03-0.1 (Rigault et al. 2025, Sect. 4.7).
Z_VOLUME_LIMIT = 0.06
Z_HUBBLE_FLOW_MIN = 0.023
Z_DR2_CLEAN = (0.03, 0.10)
Z_MIN_PLOT = 0.003

#: Error model for distance moduli: peculiar velocities of 250 km/s (the value
#: Pantheon+ assigns, VPECERR) and, for redshifts measured from the supernova's
#: own spectral features, at least 3e-3 in z (~900 km/s; Rigault et al. 2025).
SIGMA_VPEC = 250.0
SIGMA_Z_SN_FEATURES = 3e-3

# ---------------------------------------------------------------------------
# Published values, for comparison only
# ---------------------------------------------------------------------------

#: Brout et al. 2022 (ApJ 938, 110), Table 2, baseline BS21 scatter model.
PANTHEON_STANDARDISATION = {"alpha": (0.148, 0.004), "beta": (3.09, 0.04),
                            "gamma": (-0.003, 0.007), "rms": 0.171}

#: Rigault et al. 2025 (A&A 694, A1): sample sizes and Hubble-residual scatter.
DR2_PUBLISHED = {
    "n_released": 3628, "n_good_sampling": 2960, "n_basic_cuts": 2667,
    "n_volume_limited": 994, "n_non_peculiar": 2629,
    "nmad_before": 0.33, "std_before": 0.46,          # all 2629, unstandardised
    "nmad_after": 0.165, "std_after": 0.209,          # all 2629, standardised
    "nmad_after_clean": 0.150,                        # 0.03 < z < 0.1
    "std_after_clipped": 0.187, "nmad_after_clipped": 0.159, "n_after_clipped": 2606,
    "sigma_sn": 0.145,
    # what they standardised with (Ginolin et al. 2025a), local-colour step
    "alpha": 0.16, "beta": 3.05, "gamma_local_colour": 0.145,
    # selection function: P(classified | peak mag m) = 1 - S(18.8, 4.5)
    "selection_m50": 18.8, "selection_k": 4.5,
}

#: Ginolin et al. 2025a (A&A, ZTF DR2 issue; arXiv:2405.20965), Sect. 4.3:
#: 938 volume-limited (z < 0.06) SNe Ia, linear stretch, SALT2-T21; gamma for
#: the local-colour step, and (their table) for the global stellar-mass step.
#: Broken-alpha fit: two slopes either side of x1 = -0.48.
GINOLIN_2025A = {"alpha": (0.161, 0.010), "beta": (3.05, 0.06),
                 "gamma_mass": (0.145, 0.021), "gamma_local_colour": (0.143, 0.025),
                 "alpha_low": (0.271, 0.011), "alpha_high": (0.083, 0.009)}

#: Lacroix et al. 2025 (arXiv:2509.04073): scene-modelling photometry vs DR2.
#: Distances differ by a global ~0.09 mag; colours by ~0.01 mag.
LACROIX_2025_OFFSET_MAG = 0.09
