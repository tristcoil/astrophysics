"""Published numbers this project grades itself against (never fit inputs).

Every value was checked against the paper's abstract or text on 2026-10-04.
Asymmetric errors are kept as (value, minus, plus).
"""

from __future__ import annotations

#: Ghirardini et al. (2024, A&A 689, A298), eRASS1 cluster counts, flat LCDM.
GHIRARDINI2024 = {
    "Om": (0.29, 0.02, 0.01),
    "sigma8": (0.88, 0.02, 0.02),
    "S8": (0.86, 0.01, 0.01),
    "n_clusters": 5259,
    "area_deg2": 12791,
    "label": "eRASS1 clusters (Ghirardini+ 2024)",
}

#: Planck 2018 VI (Aghanim et al. 2020), TT,TE,EE+lowE+lensing, Table 2.
PLANCK2018 = {"Om": (0.3153, 0.0073, 0.0073), "sigma8": (0.8111, 0.0060, 0.0060),
              "S8": (0.832, 0.013, 0.013), "corr_Om_sigma8": None,
              "label": "Planck 2018 CMB (TT,TE,EE+lowE+lensing)"}

#: Planck 2018 best-fit point (same likelihood, Table 1), run through our CAMB
#: path as a check (published derived sigma_8 of this point: 0.8120).
PLANCK2018_BESTFIT = {"ombh2": 0.022383, "omch2": 0.12011, "theta100": 1.040909, "tau": 0.0543,
                      "As": 2.1005e-9, "ns": 0.96605, "sigma8": 0.8120}

#: Optical depth from Planck 2018 TT+lowE (Table 2), as used by project D3.
PLANCK_TAU_TT_LOWE = (0.0522, 0.0080)

#: Project D3's own fit to Planck's TT spectrum (09-planck-cmb/data/derived/
#: tt_fit.json, 2026-09-29): (ombh2, omch2, 100 theta_MC, ln 10^10 A_s, n_s),
#: tau fixed at 0.0522. Copied here so this folder runs on its own.
D3_TT_FIT = {
    "mean": [0.022210836834538488, 0.12142427668531083, 1.040683382014791, 3.0438697272087976, 0.9628705248868852],
    "err": [0.00018550774157662933, 0.002119437541954166, 0.00047493138583814653, 0.005437691990137444,
            0.0051265191629836955],
    "corr": [[1.0, -0.4852, 0.4349, -0.3893, 0.3622], [-0.4852, 1.0, -0.411, 0.9427, -0.8612],
             [0.4349, -0.411, 1.0, -0.3102, 0.4602], [-0.3893, 0.9427, -0.3102, 1.0, -0.9201],
             [0.3622, -0.8612, 0.4602, -0.9201, 1.0]],
    "tau_fixed": 0.0522,
}

#: Late-universe S8 measurements for the comparison panel.
S8_MEASUREMENTS = [
    # (label, S8, minus, plus, kind)
    ("KiDS-Legacy cosmic shear (Wright+ 2025)", 0.815, 0.021, 0.016, "lensing"),
    ("DES Y6 cosmic shear (2026)", 0.798, 0.015, 0.014, "lensing"),
    ("DES Y6 3x2pt (2026)", 0.789, 0.012, 0.012, "lensing"),
    ("SPT clusters + DES/HST lensing (Bocquet+ 2024)", 0.795, 0.029, 0.029, "clusters"),
    ("ACT DR6 CMB lensing + BAO (Madhavacheril+ 2024)", 0.840, 0.028, 0.028, "cmb-lensing"),
]

#: Omega_m of the shear surveys where quoted (DES Y6 3x2pt abstract).
DES_Y6_3X2_OM = (0.333, 0.028, 0.023)
SPT2024 = {"Om": (0.286, 0.032, 0.032), "sigma8": (0.817, 0.026, 0.026), "S8": (0.795, 0.029, 0.029)}

#: The hydrostatic mass bias: M_hydrostatic / M_true. Planck's SZ masses
#: (PSZ2 MSZ) are calibrated on X-ray hydrostatic masses with no correction,
#: so lensing-to-Planck mass ratios measure (1 - b) directly.
HYDROSTATIC_BIAS = [
    ("CMB + Planck SZ counts need (Planck 2015 XXIV)", 0.58, 0.04),
    ("Weighing the Giants lensing (von der Linden+ 2014)", 0.688, 0.072),
    ("CCCP lensing (Hoekstra+ 2015)", 0.76, 0.078),  # 0.05 stat (+) 0.06 syst
    ("Planck's default assumption (2013-2015)", 0.80, 0.0),
]
