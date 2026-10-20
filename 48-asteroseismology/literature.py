"""Published numbers we grade ourselves against (never inputs to a fit).

Star-by-star catalogue values (Yu et al. 2018, APOKASC-3, Vrard et al. 2016,
Lund et al. 2017) are read from the CDS tables by ``catalogs.py``; the
handful of numbers below come from the papers' text and tables, each with
its source.
"""

from __future__ import annotations

#: Dnu = alpha numax^beta for red giants.
DNU_NUMAX = {
    # Stello et al. 2009 (MNRAS 400, L80), their Eq. 1
    "Stello et al. 2009": {"alpha": 0.263, "alpha_err": 0.009, "beta": 0.772, "beta_err": 0.005},
    # Huber et al. 2011 (ApJ 743, 143), Table 1, SYD pipeline, numax < 300 uHz
    "Huber et al. 2011 (SYD)": {"alpha": 0.267, "alpha_err": 0.002, "beta": 0.760, "beta_err": 0.002},
}

#: The Sun's numax and Dnu as measured by other pipelines on VIRGO data.
SUN = {
    "Huber et al. 2011 (SYD)": {"numax": 3090.0, "numax_err": 30.0, "dnu": 135.1, "dnu_err": 0.1},
    # APOKASC-2 / APOKASC-3 reference values (Pinsonneault et al. 2018, 2025)
    "APOKASC-3": {"numax": 3076.0, "dnu": 135.1416},
}

#: 16 Cyg A and B.
SIXTEEN_CYG = {
    "16 Cyg A": {
        # Interferometric radius, CHARA/PAVO (White et al. 2013, MNRAS 433, 1262, Table 3)
        "radius_interf": 1.22, "radius_interf_err": 0.02,
        # Asteroseismic modelling, all constraints (Metcalfe et al. 2015, ApJ 811, L37, Table 1)
        "mass_model": 1.08, "mass_model_err": 0.02,
        "radius_model": 1.229, "radius_model_err": 0.008,
        "age_model_gyr": 7.07,
    },
    "16 Cyg B": {
        "radius_interf": 1.12, "radius_interf_err": 0.02,
        "mass_model": 1.04, "mass_model_err": 0.02,
        "radius_model": 1.116, "radius_model_err": 0.006,
        "age_model_gyr": 6.74,
    },
}

#: Bedding et al. 2011 (Nature 471, 608): observed mean period spacings of
#: the l = 1 mixed modes (Delta P_obs) for their Fig. 2 stars, and the true
#: g-mode spacing of KIC 6928997 from their Fig. 1.
BEDDING_2011 = {
    6779699: {"dp_obs": 53.0, "state": "hydrogen-shell burning (RGB)"},
    4902641: {"dp_obs": 96.0, "state": "helium-core burning"},
    6928997: {"dp_g": 77.1, "dp_obs": 55.0, "state": "hydrogen-shell burning (RGB)"},
}

#: Typical uncertainties quoted by the catalogues (for the README).
YU_TYPICAL = {"numax": 0.016, "dnu": 0.006, "mass": 0.078, "radius": 0.029}
APOKASC3_TYPICAL = {"numax": 0.006, "dnu": 0.006, "mass": 0.038, "radius": 0.018}
