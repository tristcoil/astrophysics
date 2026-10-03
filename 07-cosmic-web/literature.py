"""Published numbers we compare against -- never used to produce a result.

Everything the pipeline measures is measured from the SDSS catalogue. The
values below are consulted only afterwards, to grade those measurements, and
each carries its source.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Zehavi et al. 2011, ApJ 736, 59 -- "Galaxy Clustering in the Completed SDSS
# Redshift Survey: The Dependence on Color and Luminosity" (arXiv:1005.2413).
# SDSS DR7 main galaxy sample, NYU-VAGC, luminosity-threshold sample Mr < -20
# (0.1-band, evolved to z = 0.1, h = 1), 0.02 < z < 0.1064, 10 < r < 17.6.
# ---------------------------------------------------------------------------

ZEHAVI_REF = "Zehavi et al. 2011, ApJ 736, 59"

#: Table 2, row Mr < -20.0 (full redshift range).
ZEHAVI_MR20 = {
    "mag_limit": -20.0,
    "cz_max_kms": 31_900,
    "n_gal": 141_733,
    "n_blue": 62_862,
    "n_red": 78_868,
    "r0": 5.54, "r0_err": 0.14,       # full-covariance power-law fit, h^-1 Mpc
    "gamma": 1.83, "gamma_err": 0.01,
    "r0_diag": 6.00, "r0_diag_err": 0.09,  # diagonal-errors-only fit
    "gamma_diag": 1.79, "gamma_diag_err": 0.01,
    "number_density": 0.656e-2,        # h^3 Mpc^-3
    "pi_max": 60.0,
    "fit_range": (0.13, 20.0),         # fits use r_p < 20 h^-1 Mpc, above the 55" fibre scale
}

#: Table 2, the same threshold cut at cz < 19,250 km/s to exclude the Sloan
#: Great Wall -- which alone lowers r0 from 5.54 to 5.24.
ZEHAVI_MR20_NO_SGW = {"cz_max_kms": 19_250, "r0": 5.24, "r0_err": 0.28,
                      "gamma": 1.87, "gamma_err": 0.03}

#: Table C8: w_p(r_p) for the Mr < -20 threshold sample, h^-1 Mpc, with the
#: diagonal jackknife errors. r_p is the pair-weighted separation of each bin.
ZEHAVI_WP_RP = np.array([0.17, 0.27, 0.42, 0.67, 1.1, 1.7, 2.7, 4.2, 6.7, 10.6, 16.9, 26.8, 42.3])
ZEHAVI_WP_MR20 = np.array([366.1, 264.3, 184.0, 128.6, 84.7, 59.4, 42.9, 30.9, 21.9, 14.6,
                           8.24, 4.88, 3.58])
ZEHAVI_WP_MR20_ERR = np.array([9.3, 7.6, 6.6, 5.5, 4.3, 3.6, 3.3, 3.1, 2.7, 2.1, 1.32, 1.06, 0.85])

#: Fraction of main-sample targets left without a spectrum because two
#: fibres cannot sit closer than 55 arcsec on one plate (Zehavi et al. 2011, 2.1).
FIBRE_COLLISION_FRACTION = 0.07
FIBRE_COLLISION_ARCSEC = 55.0

# ---------------------------------------------------------------------------
# Gott et al. 2005, ApJ 624, 463 -- "A Map of the Universe" (astro-ph/0310571).
# The Sloan Great Wall, measured in a 4-degree equatorial slice, H0 = 71.
# ---------------------------------------------------------------------------

GOTT_REF = "Gott et al. 2005, ApJ 624, 463"
GOTT_SGW = {
    "ra_range_hours": (8.7, 13.7),
    "distance_range_mpc": (215.0, 370.0),
    "median_distance_mpc": 310.0,
    "redshift": 0.073,
    "comoving_length_mpc": 450.0,   # measured along the curve of the wall
    "length_mly": 1_370.0,          # the famous "1.37 billion light years" (proper, at z = 0.073)
    "h": 0.71,
}

# ---------------------------------------------------------------------------
# Strauss et al. 2002, AJ 124, 1810 -- the Main Galaxy Sample's design.
# ---------------------------------------------------------------------------

STRAUSS_REF = "Strauss et al. 2002, AJ 124, 1810"
STRAUSS_MGS = {"r_limit": 17.77, "median_z": 0.104, "surface_density_deg2": 90.0}

# ---------------------------------------------------------------------------
# Colless & Dunn 1996, ApJ 458, 435 -- "Structure and Dynamics of the Coma
# Cluster" (astro-ph/9508070). 465 members by 3-sigma clipping.
# ---------------------------------------------------------------------------

COLLESS_REF = "Colless & Dunn 1996, ApJ 458, 435"
COLLESS_COMA = {"sigma_v": 1015.0, "sigma_v_err": 59.0, "cz": 6917.0, "n_members": 465}
