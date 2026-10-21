"""Paths, lens metadata and every published number this project checks against.

Delay convention used throughout this project
---------------------------------------------
A pair (X, Y) with value d means ``d = t_Y - t_X``: image Y arrives d days
*after* image X. "RX J1131: (A, D) = +91 d" therefore reads "D lags A by 91 days".

COSMOGRAIL papers quote Delta t_XY = t_X - t_Y instead ("a negative value of
Delta t_AB means that A leads", Millon et al. 2020a), i.e. the opposite sign.
Every value below was converted from the paper's own convention, and the
paper's number is kept in a comment so the conversion can be checked by eye.
Asymmetric errors are stored as (value, err_minus, err_plus) after conversion
(a sign flip swaps the two error bars).
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIG = HERE / "figures"
RESULTS = HERE / "results"  # small JSON/CSV products (delays, H0 posteriors)

# --- Physical constants (CODATA / IAU) ---
C_KM_S = 299_792.458  # speed of light, km/s
MPC_KM = 3.085_677_581_491_367e19  # 1 Mpc in km
DAY_S = 86_400.0
ARCSEC_RAD = 4.848_136_811_095_36e-6  # 1 arcsec in radians

# ---------------------------------------------------------------------------
# Light-curve sources (all VizieR/CDS copies of the journals' tables)
# ---------------------------------------------------------------------------

#: catalogue -> (local folder, files). Files keep their CDS names because the
#: astropy CDS reader matches each data file to its ReadMe section by name.
CDS_CATALOGS = {
    # Millon et al. 2020a, A&A 640, A105 (COSMOGRAIL XIX): 15-year R-band light
    # curves of 23 lensed quasars from the 1.2 m Euler telescope; RX J1131 and
    # HE 0435 also include the Maidanak, Mercator and SMARTS data.
    "J/A+A/640/A105": ("cosmograil19", ["ReadMe", "stars.dat"]),  # + lc files, see fetch_data
    # Tewes et al. 2013b, A&A 556, A22: the 9-season RX J1131 curves behind the
    # H0LiCOW model's delays.
    "J/A+A/556/A22": ("tewes2013_rxj1131", ["ReadMe", "fig4.dat"]),
    # Bonvin et al. 2017, MNRAS 465, 4914 (H0LiCOW V): 13-season HE 0435 curves.
    "J/MNRAS/465/4914": ("bonvin2017_he0435", ["ReadMe", "fig2.dat"]),
    # Bonvin et al. 2018, A&A 616, A183 (COSMOGRAIL XVII): PG 1115+080.
    "J/A+A/616/A183": ("bonvin2018_pg1115", ["ReadMe", "fig2.dat"]),
    # Bonvin et al. 2019, A&A 629, A97 (COSMOGRAIL XVIII): WFI 2033-4723 from
    # four instruments.
    "J/A+A/629/A97": ("bonvin2019_wfi2033", ["ReadMe", "table1.dat", "c2.dat", "ecam.dat", "smarts.dat", "wfi.dat"]),
    # Eulaers et al. 2013, A&A 553, A121 (COSMOGRAIL XII): SDSS J1206+4332.
    "J/A+A/553/A121": ("eulaers2013_j1206", ["ReadMe", "j1206.dat", "hs2209.dat"]),
    # Courbin et al. 2018, A&A 609, A71 (COSMOGRAIL XVI): DES J0408-5354.
    "J/A+A/609/A71": ("courbin2018_des0408", ["ReadMe", "figure3.dat"]),
}
CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp/"

# ---------------------------------------------------------------------------
# Lens-model products
# ---------------------------------------------------------------------------

#: H0LiCOW public release (Wong et al. 2020; github.com/shsuyu/H0LiCOW-public).
H0LICOW_RAW = "https://raw.githubusercontent.com/shsuyu/H0LiCOW-public/master/"
H0LICOW_FILES = [
    "h0licow_distance_chains/RXJ1131_AO+HST_Dd_Ddt.dat",  # 11.5 MB
    "h0licow_distance_chains/HE0435_Ddt_AO+HST.dat",  # 3.6 MB
    "h0licow_distance_chains/J1206_final.csv",  # 3.4 MB
    "h0licow_distance_chains/wfi2033_dt_bic.dat",  # 20.6 MB
    "h0licow_distance_chains/PG1115_AO+HST_Dd_Ddt.dat",  # 83.6 MB
    "MontePython_cosmo_sampling/data/timedelay_6lenses/B1608_Dd_Ddt_params.dat",
    "H0_inference_code/lensutils.py",  # their likelihood code, kept for reference
    # Their own flat-LCDM chain for the six lenses, to compare posteriors with.
    "cosmo_parameter_chains/FLCDM/FLCDM_B1608DdDdt+RXJ1131_AO_HST+HE0435_AO_HST+J1206+WFI2033+PG1115.dat",
]

#: TDCOSMO 2025 public release (arXiv:2506.03023; github.com/TDCOSMO/TDCOSMO2025_public).
TDCOSMO_RAW = "https://raw.githubusercontent.com/TDCOSMO/TDCOSMO2025_public/main/"
TDCOSMO_LENS_NAMES = [
    "B1608+656", "DES0408-5354", "HE0435-1223", "PG1115+080",
    "RXJ1131-1231", "SDSS1206+4332", "WFI2033-4723", "WGD2038-4008",
]
TDCOSMO_FILES = (
    ["chains_export/" + f for f in ("LambdaCDM1a.h5", "ULambdaCDM1.h5", "LambdaCDM1d.h5", "LambdaCDM3a.h5")]
    + ["TDCOSMO_sample/%s_const_processed.pkl" % n for n in TDCOSMO_LENS_NAMES]
    + ["TDCOSMO_sample/tdcosmo_sample.yaml", "Chain_release_README.md", "likelihood_sampling.py",
       "parameter_config.yaml", "kin_axi_jam_scaling/tdcosmo_correction.pickle"]
)

# ---------------------------------------------------------------------------
# The lenses
# ---------------------------------------------------------------------------
#
# ``model_delays``: the delays that went into the published H0LiCOW/TDCOSMO lens
#   model (and hence into the public D_dt posterior). TDCOSMO 2025 (Sect. 10)
#   confirms it used "the same identical time delays and lens models" for the
#   seven lenses it shares with TDCOSMO-1/H0LiCOW.
# ``cosmograil19``: Millon et al. 2020a's re-analysis, where one exists.
# ``ddt_delay_frac``: the time-delay share of the published D_dt error budget
#   (Wong et al. 2020, Table 3).
# ``wong20_h0``: per-lens H0 in flat LCDM (Wong et al. 2020, Fig. 2 legend).

LENSES = {
    "RXJ1131": dict(
        name="RX J1131−1231", z_d=0.295, z_s=0.654,
        images=["A", "B", "C", "D"],
        lc=("cosmograil19", "lcabcd/RXJ1131_all.dat"),
        lc_h0licow=("tewes2013_rxj1131", "fig4.dat"),
        # Suyu et al. 2013, Sect. 4.2: "relative to image B: Delta t_AB = 0.7 +- 1.4,
        # Delta t_CB = -0.4 +- 2.0, Delta t_DB = 91.4 +- 1.5 days" (Tewes et al.
        # 2013b, regression difference; stat. + sys. summed linearly).
        model_delays={("B", "A"): (0.7, 1.4, 1.4), ("B", "C"): (-0.4, 2.0, 2.0), ("B", "D"): (91.4, 1.5, 1.5)},
        model_delay_ref="Tewes et al. 2013b (via Suyu et al. 2013)",
        # Millon et al. 2020a, Table 3: AB 1.6+-0.7, AC -1.0+-1.2, AD -92.5+1.9-1.8,
        # BC -2.0+1.4-1.3, BD -93.7+-2.0, CD -91.8+2.2-2.3.
        cosmograil19={("A", "B"): (-1.6, 0.7, 0.7), ("A", "C"): (1.0, 1.2, 1.2), ("A", "D"): (92.5, 1.9, 1.8),
                      ("B", "C"): (2.0, 1.4, 1.3), ("B", "D"): (93.7, 2.0, 2.0), ("C", "D"): (91.8, 2.2, 2.3)},
        h0licow="RXJ1131", tdcosmo="RXJ1131-1231",
        ddt_delay_frac=0.016, wong20_h0=(78.2, 3.4, 3.4),
    ),
    "HE0435": dict(
        name="HE 0435−1223", z_d=0.4546, z_s=1.693,
        images=["A", "B", "C", "D"],
        lc=("cosmograil19", "lcabcd/HE0435_all.dat"),
        lc_h0licow=("bonvin2017_he0435", "fig2.dat"),
        # Bonvin et al. 2017 free-knot splines, as used by Wong et al. 2017:
        # Delta t_AB = -8.8 +- 0.8, Delta t_AC = -1.1 +- 0.7, Delta t_AD = -13.8 +- 0.9.
        model_delays={("A", "B"): (8.8, 0.8, 0.8), ("A", "C"): (1.1, 0.7, 0.7), ("A", "D"): (13.8, 0.9, 0.9)},
        model_delay_ref="Bonvin et al. 2017 (free-knot splines)",
        # Millon et al. 2020a: AB -9.0+-0.8, AC -0.8+0.8-0.7, AD -13.8+-0.8,
        # BC 7.8+-0.9, BD -5.4+0.9-0.8, CD -13.2+-0.8.
        cosmograil19={("A", "B"): (9.0, 0.8, 0.8), ("A", "C"): (0.8, 0.8, 0.7), ("A", "D"): (13.8, 0.8, 0.8),
                      ("B", "C"): (-7.8, 0.9, 0.9), ("B", "D"): (5.4, 0.9, 0.8), ("C", "D"): (13.2, 0.8, 0.8)},
        h0licow="HE0435", tdcosmo="HE0435-1223",
        ddt_delay_frac=0.053, wong20_h0=(71.7, 4.5, 4.8),
        # HST/NIC2 relative astrometry (Courbin et al. 2011, Table 2), arcsec,
        # (Delta alpha, Delta delta) relative to image A; G is the lens galaxy.
        astrometry={"A": (0.0, 0.0), "B": (-1.4743, 0.5518), "C": (-2.4664, -0.6022),
                    "D": (-0.9378, -1.6160), "G": (-1.1706, -0.5665)},
        astrometry_err=0.002,  # arcsec, the systematic floor Courbin et al. add
    ),
    "PG1115": dict(
        name="PG 1115+080", z_d=0.311, z_s=1.722,
        images=["A", "B", "C"],  # A = A1 + A2 (unresolved merging pair)
        lc=("bonvin2018_pg1115", "fig2.dat"),
        # Bonvin et al. 2018 (via Millon et al. 2020a): AB -8.3+1.5-1.6,
        # AC 9.9+-1.1, BC 18.8+-1.6 (C arrives first).
        model_delays={("A", "B"): (8.3, 1.5, 1.6), ("A", "C"): (-9.9, 1.1, 1.1), ("B", "C"): (-18.8, 1.6, 1.6)},
        model_delay_ref="Bonvin et al. 2018",
        h0licow="PG1115", tdcosmo="PG1115+080",
        ddt_delay_frac=0.064, wong20_h0=(81.1, 7.1, 8.0),
    ),
    "WFI2033": dict(
        name="WFI 2033−4723", z_d=0.6575, z_s=1.662,
        images=["A1", "A2", "B", "C"],
        lc=("bonvin2019_wfi2033", ["ecam.dat", "c2.dat", "smarts.dat", "wfi.dat"]),
        # Bonvin et al. 2019: A1B 36.2+2.3-1.6, A2B 37.3+3.0-2.6, BC -59.4+1.3-3.4
        # (B arrives first).
        model_delays={("B", "A1"): (36.2, 1.6, 2.3), ("B", "A2"): (37.3, 2.6, 3.0), ("B", "C"): (59.4, 1.3, 3.4)},
        model_delay_ref="Bonvin et al. 2019",
        h0licow="WFI2033", tdcosmo="WFI2033-4723",
        ddt_delay_frac=0.029, wong20_h0=(71.6, 4.9, 3.8),
    ),
    "J1206": dict(
        name="SDSS J1206+4332", z_d=0.745, z_s=1.789,
        images=["A", "B"],
        lc=("eulaers2013_j1206", "j1206.dat"),
        # Birrer et al. 2019, Sect. 5.3: re-analysis of the Eulaers et al. 2013
        # curves, "Delta t_AB = 111.8 +2.4 -2.7 days" (Fig. 4 shows -111.8 in the
        # PyCS convention, A leading).
        model_delays={("A", "B"): (111.8, 2.4, 2.7)},
        model_delay_ref="Birrer et al. 2019 (re-analysis of Eulaers et al. 2013)",
        h0licow="J1206", tdcosmo="SDSS1206+4332",
        ddt_delay_frac=0.023, wong20_h0=(68.9, 5.1, 5.4),
    ),
    "DES0408": dict(
        name="DES J0408−5354", z_d=0.597, z_s=2.375,
        images=["A", "B", "D"],
        lc=("courbin2018_des0408", "figure3.dat"),
        # Courbin et al. 2018: AB -112.1+-2.1, AD -155.5+-12.8, BD -42.4+-17.6.
        model_delays={("A", "B"): (112.1, 2.1, 2.1), ("A", "D"): (155.5, 12.8, 12.8), ("B", "D"): (42.4, 17.6, 17.6)},
        model_delay_ref="Courbin et al. 2018",
        h0licow=None, tdcosmo="DES0408-5354",
        ddt_delay_frac=None, wong20_h0=None,
    ),
    # Not re-measured here: B1608's delays are from VLA radio monitoring
    # (Fassnacht et al. 2002) whose light curves are not public in machine-
    # readable form, and WGD 2038's curves await Dux et al. (in prep.).
    "B1608": dict(
        name="B1608+656", z_d=0.6304, z_s=1.394,
        images=["A", "B", "C", "D"],
        lc=None,
        # Fassnacht et al. 2002: B leads; AB 31.5+2.0-1.0, CB 36.0+-1.5, DB 77.0+2.0-1.0.
        model_delays={("B", "A"): (31.5, 1.0, 2.0), ("B", "C"): (36.0, 1.5, 1.5), ("B", "D"): (77.0, 1.0, 2.0)},
        model_delay_ref="Fassnacht et al. 2002 (VLA radio)",
        h0licow="B1608", tdcosmo="B1608+656",
        ddt_delay_frac=0.017, wong20_h0=(71.0, 3.3, 2.9),
    ),
    "WGD2038": dict(
        name="WGD 2038−4008", z_d=0.2283, z_s=0.777,
        images=["A", "B", "C", "D"],
        lc=None,
        model_delays=None,  # 3x3 covariance in Wong et al. 2024 (TDCOSMO XVI), Fig. 2
        model_delay_ref="Wong et al. 2024",
        h0licow=None, tdcosmo="WGD2038-4008",
        ddt_delay_frac=None, wong20_h0=None,
    ),
}

#: The lenses whose delays we re-measure ourselves, in display order.
MEASURED = ["RXJ1131", "HE0435", "PG1115", "WFI2033", "J1206", "DES0408"]
#: The six H0LiCOW lenses (Wong et al. 2020).
H0LICOW_SIX = ["B1608", "RXJ1131", "HE0435", "J1206", "WFI2033", "PG1115"]

# ---------------------------------------------------------------------------
# Published H0 values we compare with (km/s/Mpc; value, err_minus, err_plus)
# ---------------------------------------------------------------------------

H0_LIT = {
    # Wong et al. 2020 (H0LiCOW XIII), flat LCDM, six lenses, Table 6.
    "H0LiCOW 2020": (73.3, 1.8, 1.7),
    # Millon et al. 2020b (TDCOSMO I), seven lenses, power-law models.
    "TDCOSMO-I 2020": (74.2, 1.6, 1.6),
    # Birrer et al. 2020 (TDCOSMO IV), seven lenses, mass-sheet free (kinematics only).
    "TDCOSMO-IV 2020": (74.5, 6.1, 5.6),
    # TDCOSMO Collaboration 2025 (arXiv:2506.03023), Table 6.
    "TDCOSMO 2025": (71.6, 3.3, 3.9),  # LCDM1a: 8 lenses + Pantheon+ Omega_m
    "TDCOSMO 2025 + SLACS + SL2S": (74.3, 3.7, 3.1),  # LCDM1d
    "TDCOSMO 2025 alone": (73.7, 4.4, 4.7),  # ULCDM1, uniform Omega_m
}

#: Other routes to H0 drawn on the comparison figure. D1 and D3 are this
#: repository's own measurements; C4 and I9 are not built yet, so their
#: published values are used and labelled as such.
H0_OTHER = {
    "D1 Cepheids + SNe Ia (this repo)": (73.55, 1.02, 1.02),
    "D3 Planck CMB (this repo)": (66.66, 0.90, 0.90),
    "C4 GW170817 standard siren (Abbott+2017)": (70.0, 8.0, 12.0),
    "I9 GWTC-5.0 dark + bright sirens (LVK 2026)": (71.0, 7.1, 9.0),
}
