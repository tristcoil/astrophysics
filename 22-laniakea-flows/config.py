"""Project configuration: paths, data sources, reference frames, and every assumption.

Every number this project *chooses* rather than measures lives here: the
reference frames that turn heliocentric redshifts into cosmological ones, the
cosmology that turns redshifts into distances, the Hubble constant that turns
distances into peculiar velocities, the radial bins of the bulk-flow fits, the
smoothing of the density field. Each sits next to the reason for it and, where
it comes from a paper, the paper. Published results we compare against are
collected at the bottom, so the analysis code never hides a literature value.

Distances are in h^-1 Mpc (H0 = 100 h km/s/Mpc) wherever the Hubble constant
should not matter -- maps, bulk flows, the density field -- so that 1 h^-1 Mpc
is 100 km/s of Hubble expansion. Only the Hubble diagram uses physical Mpc.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
CF4_DIR = DATA_DIR / "cf4"
TWOMRS_DIR = DATA_DIR / "2mrs"
GROUPS_DIR = DATA_DIR / "tully2015"
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, CF4_DIR, TWOMRS_DIR, GROUPS_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-project-atlas/J4 (hobby research; python-urllib)"
C_KMS = 299_792.458


# --------------------------------------------------------------------------
# Data sources
# --------------------------------------------------------------------------

#: CDS FTP. Its bot challenge targets browser user agents; a plain one gets
#: the files directly.
CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp"


@dataclass(frozen=True)
class Source:
    """One file we download, and why."""

    local: str  # path under data/, as the analysis reads it (gz files are unpacked)
    url: str
    what: str
    reference: str
    gzipped: bool = False


def _cds(catalog: str, filename: str, folder: str, what: str, reference: str) -> Source:
    gz = filename.endswith(".gz")
    local = f"{folder}/{filename[:-3] if gz else filename}"
    return Source(local, f"{CDS_FTP}/{catalog}/{filename}", what, reference, gzipped=gz)


CF4_REF = "Tully et al. 2023, ApJ 944, 94 (Cosmicflows-4), VizieR J/ApJ/944/94"
TWOMRS_REF = "Huchra et al. 2012, ApJS 199, 26 (2MASS Redshift Survey), VizieR J/ApJS/199/26"
GROUPS_REF = "Tully 2015, AJ 149, 171 (2MASS galaxy groups), VizieR J/AJ/149/171"
HOWLETT_REF = "Howlett et al. 2022, MNRAS, doi:10.1093/mnras/stac1681 (SDSS peculiar velocity catalogue), Zenodo record 6824749 v1.1.0"
ZENODO_SDSS = "https://zenodo.org/api/records/6824749/files"

SOURCES: list[Source] = [
    _cds("J/ApJ/944/94", "ReadMe", "cf4", "byte-by-byte description of the three CF4 tables", CF4_REF),
    _cds("J/ApJ/944/94", "table2.dat.gz", "cf4",
         "distances of 55,877 individual galaxies, combined and per method (TF, FP, SN Ia, SBF, SN II, "
         "TRGB, Cepheids, masers)", CF4_REF),
    _cds("J/ApJ/944/94", "table3.dat.gz", "cf4",
         "distances of 38,053 groups, combined and per method, with the calibrator zero point", CF4_REF),
    _cds("J/ApJ/944/94", "table4.dat.gz", "cf4",
         "group distances, velocities in three frames, and CF4's own peculiar velocities", CF4_REF),
    _cds("J/ApJS/199/26", "ReadMe", "2mrs", "byte-by-byte description of the 2MRS tables", TWOMRS_REF),
    _cds("J/ApJS/199/26", "table3.dat.gz", "2mrs",
         "the 44,599 galaxies of the 2MASS Redshift Survey (Ks <= 11.75): positions, magnitudes, redshifts",
         TWOMRS_REF),
    _cds("J/AJ/149/171", "ReadMe", "tully2015", "byte-by-byte description of the group tables", GROUPS_REF),
    _cds("J/AJ/149/171", "table5.dat.gz", "tully2015",
         "the 43,038 2MRS galaxies with their group ('nest') memberships, used to collapse fingers of God",
         GROUPS_REF),
    # Validation only: the SDSS Fundamental Plane catalogue that CF4's largest
    # component comes from, to check how its distances entered CF4.
    Source("howlett2022/SDSS_PV_public.dat", f"{ZENODO_SDSS}/SDSS_PV_public.dat/content",
           "34,059 SDSS early-type galaxies with Fundamental Plane log-distance ratios", HOWLETT_REF),
    Source("howlett2022/data_description.pdf", f"{ZENODO_SDSS}/data_description.pdf/content",
           "column description of the SDSS PV catalogue", HOWLETT_REF),
]


# --------------------------------------------------------------------------
# Reference frames
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Motion:
    """A velocity in Galactic coordinates, with its published uncertainties."""

    name: str
    speed: float  # km/s
    l: float  # deg
    b: float  # deg
    e_speed: float
    e_l: float
    e_b: float
    reference: str


#: The Sun relative to the cosmic microwave background: the Solar dipole of
#: Planck 2018 (HFI), Planck Collaboration 2020, A&A 641, A1, Tables 2 and 3:
#: 3362.08 +- 0.99 uK, i.e. v/c = (1.23357 +- 0.00036) x 10^-3.
SUN_CMB = Motion("Sun relative to the CMB", 369.82, 264.021, 48.253, 0.11, 0.011, 0.005,
                 "Planck Collaboration 2020, A&A 641, A1, Table 3 (HFI 2018 dipole)")

#: The Sun relative to the Local Group's barycentre (Diaz et al. 2014), as
#: adopted in the same Planck table.
SUN_LG = Motion("Sun relative to the Local Group", 299.0, 98.4, -5.9, 15.0, 3.6, 3.0,
                "Diaz et al. 2014, via Planck Collaboration 2020, A&A 641, A1, Table 3")

#: The Local Group relative to the CMB: the vector difference of the two above,
#: as tabulated by Planck (620 +- 15 km/s). Our own vector subtraction is
#: checked against it in ``catalogs.local_group_motion``.
LG_CMB = Motion("Local Group relative to the CMB", 620.0, 271.9, 29.6, 15.0, 2.0, 1.4,
                "Planck Collaboration 2020, A&A 641, A1, Table 3")


# --------------------------------------------------------------------------
# Cosmology
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Cosmology:
    """Flat Lambda-CDM with everything the analysis needs."""

    name: str
    H0: float  # km/s/Mpc
    Om0: float
    Ob0h2: float
    Oc0h2: float
    ns: float
    sigma8: float
    tau: float
    mnu: float  # eV, one massive neutrino
    reference: str

    @property
    def h(self) -> float:
        return self.H0 / 100.0


#: Planck 2018 TT,TE,EE+lowE+lensing+BAO (Planck Collaboration 2020, A&A 641,
#: A6, Table 2, last column). Nusser & Tully (2026) adopt exactly these numbers
#: (H0 67.66, Omega_m 0.3111, sigma8 0.8102, ns 0.9665) for the Lambda-CDM
#: prediction of bulk flows, so our comparison uses them too. Only the *shape*
#: of the expansion history (Omega_m) enters distances in h^-1 Mpc.
PLANCK18 = Cosmology("Planck 2018", 67.66, 0.3111, 0.02242, 0.11933, 0.9665, 0.8102, 0.0561, 0.06,
                     "Planck Collaboration 2020, A&A 641, A6, Table 2 (TT,TE,EE+lowE+lensing+BAO)")

#: Cosmicflows-4's convention for converting distances into peculiar velocities
#: (Tully et al. 2023, Sect. 13 and Table 4): H0 = 75 km/s/Mpc, and a
#: cosmological curvature correction f(z) = 1 + (1 - q0) z / 2 - (1 - q0 - 3 q0^2 + j0) z^2 / 6
#: with Omega_m = 0.27, Omega_Lambda = 0.73, so q0 = Om/2 - OL = -0.595, j0 = 1.
CF4_H0_VPEC = 75.0
CF4_OMEGA_M = 0.27
CF4_Q0 = CF4_OMEGA_M / 2 - (1 - CF4_OMEGA_M)
CF4_J0 = 1.0
#: CF4 switches from the Davis & Scrimgeour estimator nearby to Watkins &
#: Feldman's beyond this Local-Sheet velocity, with a linear ramp (their Eq. 11).
CF4_RAMP_KMS = 3000.0


# --------------------------------------------------------------------------
# 1. The Hubble diagram
# --------------------------------------------------------------------------

#: Tully et al. (2023) average log10 H_i over groups beyond 4,000 km/s, where
#: peculiar velocities (~300 km/s) are a small fraction of the expansion.
HUBBLE_VMIN = 4000.0
#: Velocity bins for the Hubble-parameter-versus-distance plot (their Fig. 18
#: uses 1,000 km/s bins).
HUBBLE_BIN_KMS = 1000.0
#: "Well-measured" groups: distance-modulus error at most 0.10 mag (5% in
#: distance), CF4's highlighted subset.
HUBBLE_GOOD_EDM = 0.10


# --------------------------------------------------------------------------
# 2. Peculiar-velocity maps
# --------------------------------------------------------------------------

#: Shells for the sky maps, in h^-1 Mpc of redshift distance (CMB frame).
#: 0-30 shows the Virgo infall and the Local Void; 30-60 the Great Attractor
#: and Perseus-Pisces; 60-100 Coma and the Hercules region; 100-160 Shapley.
SKY_SHELLS = [(0.0, 30.0), (30.0, 60.0), (60.0, 100.0), (100.0, 160.0)]
#: Smoothing of the sky maps: each HEALPix pixel shows the error-weighted mean
#: peculiar velocity of the groups within this angle of it.
SKY_SMOOTH_DEG = 12.0
SKY_NSIDE = 32
#: The supergalactic-plane slab: |SGZ| below this, in h^-1 Mpc.
SG_SLAB_HALF = 15.0


# --------------------------------------------------------------------------
# 3. Bulk flows (the method of Nusser & Tully 2026)
# --------------------------------------------------------------------------

#: Their radial bins in redshift distance, h^-1 Mpc (their Table 1). Objects
#: below cz_CMB = 4,000 km/s are never fitted: there, peculiar velocities are
#: no longer small against the expansion.
BULK_BINS = [(40.0, 60.0), (60.0, 90.0), (90.0, 120.0), (120.0, 160.0), (160.0, 200.0), (200.0, 300.0)]
BULK_CZ_MIN = 4000.0
#: A velocity floor added in quadrature to every object's error, so that the
#: few objects with tiny distance errors cannot dominate (their Eq. 7).
BULK_SIGMA_FLOOR = 100.0
#: A bin is used only if its effective number of equally weighted objects,
#: (sum w)^2 / sum w^2, is at least this (their Eq. 2).
BULK_NEFF_MIN = 10.0
#: Thin shells for "bulk flow in spheres of growing radius": the cumulative
#: sphere of radius R combines every shell inside R.
SPHERE_SHELL = 10.0
SPHERE_RMAX = 300.0
#: Linear power spectrum and velocity correlation functions (their Appendix B):
#: k from 1e-4 to 10 h/Mpc, correlation functions tabulated to 700 h^-1 Mpc
#: (the largest separation two objects within 300 h^-1 Mpc can have is 600).
PK_KMIN, PK_KMAX, PK_NK = 1e-4, 10.0, 4096
PSI_SMAX, PSI_NS = 700.0, 7001


# --------------------------------------------------------------------------
# 4. The 2MRS density field, its gravity, and the reconstruction
# --------------------------------------------------------------------------

#: 2MRS: Ks(isophotal, extinction corrected) <= 11.75, |b| >= 5 deg, and
#: |b| >= 8 deg toward the bulge (|l| < 30 deg) (Huchra et al. 2012).
TWOMRS_KLIM = 11.75
ZOA_B = 5.0
ZOA_B_BULGE = 8.0
ZOA_L_BULGE = 30.0
#: Galaxies nearer than this (in the Local Group frame) are Local Group members
#: or so close that their redshift says nothing about their distance.
TWOMRS_CZ_MIN = 300.0

#: K-band k-correction, k(z) = -6 log10(1 + z), independent of galaxy type
#: (Kochanek et al. 2001, used by Erdogdu et al. 2006 and Lavaux & Hudson 2011).
KCORR_COEFF = -6.0

#: Luminosity dependence of galaxy bias measured from 2MASS: b/b* = 0.73 + 0.24 L/L*
#: (Westover 2007, as used by Carrick et al. 2015). Brighter galaxies, which
#: dominate the far parts of a flux-limited survey, cluster more strongly.
BIAS_A, BIAS_B = 0.73, 0.24

#: The reconstruction box: a cube of side 2 x RECON_HALF centred on us, with
#: galaxies used out to RECON_RMAX (beyond it the field is set to the mean).
#: 2MRS is too sparse beyond ~200 h^-1 Mpc: a galaxy there must be 4x brighter
#: than L* to make the catalogue.
RECON_HALF = 200.0
RECON_NGRID = 128  # 3.125 h^-1 Mpc cells
RECON_RMAX = 160.0
#: Gaussian smoothing for linear theory: 4 h^-1 Mpc, the scale Carrick et al.
#: (2015) found best reproduces simulated velocities.
RECON_SMOOTH = 4.0
#: Gravity is switched on gradually (Yahil et al. 1991; Carrick et al. 2015):
#: beta* rises in these steps, each iteration starting from the last positions.
RECON_BETAS = [round(0.02 * i, 2) for i in range(0, 41)]  # 0 ... 0.80
#: Positions are averaged over the last few iterations to damp oscillations
#: in triple-valued regions near clusters (Carrick et al. 2015).
RECON_AVERAGE_LAST = 5

#: The velocity-velocity comparison: CF4 groups out to this redshift distance,
#: where the 2MRS field is well sampled on all sides.
VCOMP_RMAX = 120.0
VCOMP_CZ_MIN = 300.0

#: Basins of attraction: the gravitational potential is smoothed further, so
#: that shot noise does not split basins into fragments. Results at several
#: scales show how the answer depends on it.
BASIN_SMOOTH = [6.0, 10.0, 15.0]
BASIN_SMOOTH_MAIN = 10.0


# --------------------------------------------------------------------------
# Landmarks, for labelling maps (positions J2000; heliocentric cz)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Landmark:
    name: str
    ra: float
    dec: float
    cz_helio: float  # km/s; mean for clusters
    note: str


LANDMARKS = [
    Landmark("Virgo", 187.706, 12.391, 1088.0, "M87; Virgo cluster mean velocity (Mei et al. 2007 distance 16.5 Mpc)"),
    Landmark("Fornax", 54.621, -35.451, 1425.0, "NGC 1399"),
    Landmark("Centaurus", 192.205, -41.311, 3418.0, "Abell 3526 (NGC 4696)"),
    Landmark("Hydra", 159.178, -27.528, 3777.0, "Abell 1060 (NGC 3311)"),
    Landmark("Norma", 243.887, -60.908, 4871.0, "Abell 3627, the Great Attractor's core (Woudt et al. 2008)"),
    Landmark("Perseus", 49.951, 41.512, 5366.0, "Abell 426 (NGC 1275)"),
    Landmark("Coma", 194.953, 27.981, 6925.0, "Abell 1656"),
    Landmark("Ophiuchus", 258.116, -23.363, 8400.0, "Ophiuchus cluster, z = 0.028"),
    Landmark("Hercules", 241.31, 17.72, 10970.0, "Abell 2151, z = 0.0366"),
    Landmark("Shapley", 201.987, -31.496, 14390.0, "Abell 3558, the Shapley Concentration's core, z = 0.048"),
]


# --------------------------------------------------------------------------
# Published results we compare against (never used as inputs)
# --------------------------------------------------------------------------

PUBLISHED = {
    "cf4_H0": {"value": 74.6, "stat": 0.8, "sys": 3.0, "log_H0": 1.8725, "log_rms": 0.091,
               "ref": "Tully et al. 2023, ApJ 944, 94, Sect. 12 (groups beyond 4,000 km/s)"},
    "cf4_vpec_sigma": {"value": 300.0, "ref": "Tully et al. 2023: groups with e(DM) <= 0.10 mag"},
    "watkins2023": {"V": 395.0, "err": 29.0, "R": 150.0, "p": "< 0.03%",
                    "ref": "Watkins et al. 2023, MNRAS, doi:10.1093/mnras/stad1984 (minimum variance, 150 h^-1 Mpc)"},
    "whitford2023_mle": {"V": 408.0, "err": 165.0, "R": 49.0,
                         "ref": "Whitford, Howlett & Davis 2023, MNRAS, doi:10.1093/mnras/stad2764 (MLE)"},
    "whitford2023_mve": {"V": 428.0, "err": 108.0, "R": 173.0, "p": "0.11%",
                         "ref": "Whitford, Howlett & Davis 2023, MNRAS, doi:10.1093/mnras/stad2764 (MVE)"},
    "nusser_tully2026": {
        "ref": "Nusser & Tully 2026, arXiv:2608.14265",
        "bins_nobj": {"0-40": 4227, "40-60": 4239, "60-90": 6585, "90-120": 6767, "120-160": 8481,
                      "160-200": 6347, "200-300": 19101},
        "bin_120_160": {"Vdip": 628.0, "err": 82.0, "V_SGX": -556.0, "pV": 6.67e-4, "pmV": 1.11e-4},
        "per_bin_pmV": {"40-60": 0.753, "60-90": 0.275, "90-120": 0.00209, "120-160": 1.11e-4},
        "per_bin_pV": {"40-60": 0.646, "60-90": 0.197, "90-120": 0.00208, "120-160": 6.67e-4},
        "per_bin_pm": {"40-60": 0.630, "60-90": 0.466, "90-120": 0.148, "120-160": 0.0148},
        # Figure 5: standardized deviations (m, SGX, SGY, SGZ) for All Individual
        "per_bin_sigma": {"40-60": (0.48, -1.11, -0.01, -0.67), "60-90": (0.73, -2.01, -0.20, -0.78),
                          "90-120": (1.45, -3.58, -0.79, -1.04), "120-160": (2.44, -3.53, -1.64, -1.33),
                          "160-200": (2.36, -1.92, -0.42, -1.12), "200-300": (2.09, -2.22, -1.40, -0.64)},
        "global_40_160": {"nobj": 26072, "pm": 0.185, "pV": 0.0195, "pmV": 0.0156},
        "global_40_300": {"nobj": 51520, "pm": 0.143, "pV": 0.0222, "pmV": 0.0181},
    },
    "courtois2025": {"V": 315.0, "err": 40.0, "R": 150.0,
                     "ref": "Courtois et al. 2025, A&A 701, A187 (CF4++ field)"},
    "erdogdu2006": {"misalign_50": (12.0, 7.0), "misalign_130": (21.0, 8.0), "beta_L": (0.40, 0.09),
                    "rho_L": (7.67e8, 1.02e8),
                    "ref": "Erdogdu et al. 2006, MNRAS 368, 1515 (2MRS flux dipole, Ks < 11.25)"},
    "carrick2015": {"beta": (0.431, 0.021), "Vext": (159.0, 23.0), "Vext_l": (304.0, 11.0),
                    "Vext_b": (6.0, 13.0), "fsigma8": (0.427, 0.027),
                    "ref": "Carrick et al. 2015, MNRAS 450, 317 (2M++ vs SFI++/SNe)"},
    "watkins2025": {"beta": (0.31, 0.01), "H0_vpec": (75.0, 0.1),
                    "ref": "Watkins, Clark & Feldman 2025, arXiv:2512.03168"},
    "laniakea": {"diameter_kms": 12000.0, "diameter_mpc": 160.0, "mass": 1e17,
                 "ref": "Tully, Courtois, Hoffman & Pomarede 2014, Nature 513, 71"},
    "valade2024": {"ref": "Valade et al. 2024, Nature Astronomy, doi:10.1038/s41550-024-02370-0",
                   "finding": "slight probabilistic preference for Laniakea being part of the Shapley basin"},
    "tully2008_lg": {"virgo": 185.0, "local_void": 259.0, "residual": 455.0, "total_ls": 631.0,
                     "ref": "Tully et al. 2008, ApJ 676, 184 (Local Sheet frame)"},
}
