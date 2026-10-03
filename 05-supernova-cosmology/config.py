"""Project configuration: paths, data sources, sample cuts, and reference values.

Everything tunable lives here, together with the published numbers we grade our
own results against. Those published values are used only for comparison at the
very end -- never as inputs to a fit.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
CHAIN_DIR = DATA_DIR / "chains"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, CHAIN_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: The Pantheon+ / SH0ES public data release on GitHub. Pantheon+ is the
#: compilation of 1701 Type Ia supernova light curves described in Scolnic et al.
#: (2022) and analysed for cosmology in Brout et al. (2022).
REPO = "PantheonPlusSH0ES/DataRelease"
REPO_RAW = f"https://raw.githubusercontent.com/{REPO}/main/"

_DIST = "Pantheon+_Data/4_DISTANCES_AND_COVAR/"
_CHAINS = "Pantheon+_Data/5_COSMOLOGY/chains/"

#: The data we fit. Local name -> path inside the repository.
DATA_FILES = {
    # One row per light curve: redshifts, standardised magnitudes, Cepheid host
    # distances for the calibrators, SALT2 fit parameters.
    "Pantheon+SH0ES.dat": _DIST + "Pantheon+SH0ES.dat",
    # 1701 x 1701 covariance of the standardised magnitudes: statistical errors
    # plus every systematic the team modelled (calibration, dust, peculiar
    # velocities, ...). The authors are explicit that fits must use this.
    "Pantheon+SH0ES_STAT+SYS.cov": _DIST + "Pantheon+SH0ES_STAT+SYS.cov",
    # Statistical part only. Used just to show how much systematics matter.
    "Pantheon+SH0ES_STATONLY.cov": _DIST + "Pantheon+SH0ES_STATONLY.cov",
    # Column definitions and the survey-ID key, kept alongside for provenance.
    "README_distances.txt": _DIST + "README",
}

#: The Pantheon+ team's own MCMC chains. Not used in any fit: they are the
#: answer key our independent likelihood is checked against.
CHAIN_FILES = {
    "flat_lcdm": _CHAINS + "Pantheon+Only/Pantheon+only_FlatLambdaCDM.txt.gz",
    "lcdm": _CHAINS + "Pantheon+Only/Pantheon+only_LambdaCDM.txt.gz",
    "flat_wcdm": _CHAINS + "Pantheon+Only/Pantheon+only_FlatwCDM.txt.gz",
    "shoes_flat_lcdm": _CHAINS + "Pantheon+SH0ES/Pantheon+SH0ES_FlatLambdaCDM.txt.gz",
}

#: Survey key, copied from the data release README (``IDSURVEY`` column).
SURVEYS = {
    1: "SDSS", 4: "SNLS", 5: "CSP", 10: "DES", 15: "PS1MD", 18: "CNIa0.02",
    50: "LOWZ/JRK07", 51: "LOSS1", 56: "SOUSA", 57: "LOSS2", 61: "CFA1",
    62: "CFA2", 63: "CFA3S", 64: "CFA3K", 65: "CFA4p2", 66: "CFA4p3",
    100: "HST", 101: "SNAP", 106: "CANDELS", 150: "FOUND",
}

#: Surveys grouped by how they found their supernovae, which is also roughly
#: which redshift range they cover. ``SNAP`` (101) is a label for six z = 1.0-1.4
#: supernovae observed with HST, not the satellite of that name, which never flew.
SURVEY_GROUPS = {
    "Nearby searches (CfA, CSP, LOSS, Foundation, ...)":
        (5, 18, 50, 51, 56, 57, 61, 62, 63, 64, 65, 66, 150),
    "Deep ground-based surveys (SDSS, PS1, SNLS, DES)": (1, 4, 10, 15),
    "Hubble Space Telescope": (100, 101, 106),
}

# ---------------------------------------------------------------------------
# Sample cuts -- identical to the published analysis
# ---------------------------------------------------------------------------

#: Supernovae nearer than this are dropped from cosmology fits. Galaxies have
#: their own "peculiar" motions of a few hundred km/s on top of the cosmic
#: expansion; at z = 0.01 (cz ~ 3000 km/s) that is already a ~10% error on the
#: velocity, and it gets worse closer in. Brout et al. (2022) use z > 0.01.
Z_MIN_COSMOLOGY = 0.01

#: Hubble-flow window used by SH0ES (Riess et al. 2022) for its H0 measurement.
#: Far enough that peculiar velocities are small, near enough that the result
#: barely depends on the cosmological model. Used for a cross-check only.
Z_SH0ES_HUBBLE_FLOW = (0.023, 0.15)

#: Speed of light, km/s (exact, by definition of the metre).
C_KM_S = 299_792.458

#: H0 used when it cannot be measured (fits without the Cepheid calibration).
#: With supernovae alone H0 is perfectly degenerate with the supernova absolute
#: magnitude, so any value gives identical cosmology results.
H0_FIDUCIAL = 70.0

# ---------------------------------------------------------------------------
# Parameter grids for the likelihood scans
# ---------------------------------------------------------------------------

#: Flat LambdaCDM: one free shape parameter. Same range as the published chains'
#: flat prior, 0.1 < Omega_m < 0.9.
GRID_FLAT_OM = (0.1, 0.9, 1601)

#: Non-flat LambdaCDM, the classic (Omega_m, Omega_Lambda) plane. Wide enough to
#: show Einstein-de Sitter (1, 0) and the "no Big Bang" corner of the plane.
GRID_OM = (0.0, 1.2, 241)
GRID_OL = (-0.4, 1.6, 401)

#: Coarser version of the same plane, for the covariance comparison (stat-only,
#: diagonal-only), where only the size of the contours matters.
GRID_OM_COARSE = (0.0, 1.2, 121)
GRID_OL_COARSE = (-0.4, 1.6, 201)

#: Flat wCDM: dark energy with a constant equation of state w = P / (rho c^2).
#: The grid covers where the likelihood lives; marginal posteriors are then
#: computed inside the published chains' prior box (below), because along the
#: Omega_m-w degeneracy the answer depends on where the prior cuts it off.
GRID_W_OM = (0.0, 0.7, 141)
GRID_W_W = (-2.5, -0.3, 221)

#: Flat priors used by the Pantheon+ chains (read from their headers). Our
#: marginal posteriors apply the same box so the comparison is like-for-like.
CHAIN_PRIOR_OM = (0.1, 0.9)
CHAIN_PRIOR_W = (-2.0, 0.5)

# ---------------------------------------------------------------------------
# Published values (for comparison only)
# ---------------------------------------------------------------------------

#: Brout et al. 2022, ApJ 938, 110, Table 3. Pantheon+ supernovae alone unless
#: marked SH0ES. Asymmetric errors are given as (value, +err, -err).
PUBLISHED_PANTHEON = {
    "flat_lcdm_Om": (0.334, 0.018, 0.018),
    "lcdm_Om": (0.306, 0.057, 0.057),
    "lcdm_OL": (0.625, 0.084, 0.084),
    "flat_wcdm_Om": (0.309, 0.063, 0.069),
    "flat_wcdm_w": (-0.90, 0.14, 0.14),
    "shoes_flat_lcdm_H0": (73.6, 1.1, 1.1),
    "shoes_flat_lcdm_Om": (0.334, 0.018, 0.018),
}

#: Riess et al. 2022, ApJL 934, L7 -- the SH0ES Cepheid distance ladder.
SH0ES_R22_H0 = (73.04, 1.04)

#: Planck Collaboration 2020, A&A 641, A6: base LambdaCDM from the cosmic
#: microwave background (TT,TE,EE+lowE+lensing). An *early-universe* inference
#: of H0 -- a prediction from the physics of the young universe, not a
#: measurement of the local expansion.
PLANCK18_H0 = (67.36, 0.54)
PLANCK18_OM = (0.3153, 0.0073)

#: The 1998 discovery papers, for the history lesson.
#: Riess et al. 1998 (AJ 116, 1009): Omega_Lambda > 0 at 3.0 sigma and 4.0 sigma
#: (two light-curve fitting methods), from 16 high-z + 34 nearby supernovae.
#: Perlmutter et al. 1999 (ApJ 517, 565): P(Lambda > 0) = 99% from 42 high-z
#: + 18 nearby supernovae.
DISCOVERY_1998 = {
    "riess98_sigma_OL_positive": (3.0, 4.0),
    "perlmutter99_prob_L_positive": 0.99,
}
