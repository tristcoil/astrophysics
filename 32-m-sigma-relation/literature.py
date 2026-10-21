"""Published numbers we grade ourselves against. Never used as fit inputs.

Relations are stored as log M = a + b x with x = log(sigma / 200 km/s) or
log(M_bulge / 1e11 Msun), and eps the intrinsic scatter in log M (dex).
"""

from __future__ import annotations

import math

#: Kormendy & Ho (2013, ARA&A 51, 511), Section 6.6. Their zero points are for
#: log(M / 1e9 Msun); +9 converts them to log M.
KH13 = {
    # Eq. 3: symmetric least squares (Tremaine et al. 2002), every galaxy given
    # the same log M error of 0.117 dex. Their adopted relation.
    "msigma_eq3": {"a": 9 - 0.510, "ea": 0.049, "b": 4.377, "eb": 0.290, "eps": 0.29, "eeps": 0.03,
                   "label": "Kormendy & Ho 2013 (eq. 3)"},
    # Eq. 5: the same with individual errors.
    "msigma_eq5": {"a": 9 - 0.501, "ea": 0.049, "b": 4.414, "eb": 0.295, "eps": 0.28},
    # Eq. 12: only the 29 galaxies with sigma_e < 270 km/s.
    "msigma_eq12": {"a": 9 - 0.54, "ea": 0.07, "b": 4.26, "eb": 0.44, "eps": 0.30, "n": 29},
    # Eq. 10: M / 1e9 = 0.49 (+0.06 / -0.05) (M_bulge / 1e11)^(1.16 +/- 0.08), 0.29 dex.
    "mbulge_eq10": {"a": 9 + math.log10(0.49), "ea": 0.05, "b": 1.16, "eb": 0.08, "eps": 0.29},
    "ratio_percent": {"value": 0.49, "hi": 0.06, "lo": 0.05},
}

#: McConnell & Ma (2013, ApJ 764, 184), Table 2 (MPFITEXY unless noted).
MM13 = {
    "all": {"a": 8.32, "ea": 0.05, "b": 5.64, "eb": 0.32, "eps": 0.38, "n": 72,
            "label": "McConnell & Ma 2013"},
    "all_linmix": {"a": 8.31, "ea": 0.06, "b": 5.67, "eb": 0.33, "eps": 0.40, "eeps": 0.04, "n": 72},
    "with_upper_limits_linmix": {"a": 8.15, "ea": 0.05, "b": 5.58, "eb": 0.30, "eps": 0.43, "n": 164},
    "early": {"a": 8.39, "ea": 0.06, "b": 5.20, "eb": 0.36, "eps": 0.34, "n": 53},
    "late": {"a": 8.07, "ea": 0.21, "b": 5.06, "eb": 1.16, "eps": 0.46, "n": 19},
    "mbulge_dynamical": {"a": 8.46, "ea": 0.08, "b": 1.05, "eb": 0.11, "eps": 0.34, "n": 35},
}

#: van den Bosch (2016, ApJ 831, 134), eq. 1: LINMIX_ERR with upper limits.
VDB16 = {
    "msigma": {"a": 8.32, "ea": 0.04, "b": 5.35, "eb": 0.23, "eps": 0.49, "eeps": 0.03, "n": 230,
               "label": "van den Bosch 2016"},
}

#: Virial factors <f> offered by the AGN Black Hole Mass Database (Bentz & Katz
#: 2015), each calibrated by putting reverberation-mapped AGN on a quiescent
#: M-sigma relation (Graham et al. 2011 used a relation with a lower zero point).
F_FACTOR = {
    "Onken+04": (5.5, 1.8),
    "Park+12": (5.1, 1.3),
    "Batiste+17": (4.8, 1.7),
    "Grier+13": (4.3, 1.1),
    "Graham+11": (2.8, 0.6),
}

#: Ho & Kim (2014, ApJ 789, 17): <f> calibrated separately for active galaxies
#: whose bulges are classical (or elliptical) and for pseudobulges.
HO_KIM_2014 = {"classical": (6.3, 1.5), "pseudo": (3.2, 0.7)}

#: The mass of M87's black hole, by method. All verified against the papers;
#: distances differ, and dynamical masses scale in proportion to distance.
M87_HISTORY = [
    {"ref": "Gebhardt et al. 2011", "method": "stellar dynamics", "M": 6.6e9, "hi": 0.4e9, "lo": 0.4e9,
     "D": 17.9},
    {"ref": "Walsh et al. 2013", "method": "gas dynamics", "M": 3.5e9, "hi": 0.9e9, "lo": 0.7e9, "D": 17.9},
    {"ref": "EHT Collaboration 2019", "method": "shadow (EHT)", "M": 6.5e9, "hi": 0.7e9, "lo": 0.7e9,
     "D": 16.8},
    {"ref": "Liepold et al. 2023", "method": "stellar dynamics (triaxial)", "M": 5.37e9, "hi": 0.37e9,
     "lo": 0.25e9, "D": 16.8},
]

#: Sgr A*: GRAVITY Collaboration (2022, A&A 657, L12).
GRAVITY_2022 = {"M": 4.297e6, "eM": 0.012e6, "R0": 8.277}

#: Shankar et al. (2016, MNRAS 460, 3119): most black hole hosts have higher
#: sigma than SDSS early-type galaxies of the same stellar mass, which they
#: attribute to the need to resolve the sphere of influence. Kormendy (2020,
#: IAU Symp. 353, arXiv:1909.10821) argues the relations for ellipticals and
#: classical bulges are not biased. Qualitative; we test it with ATLAS3D.
SELECTION_DEBATE = ("Shankar et al. 2016", "Kormendy 2020")
