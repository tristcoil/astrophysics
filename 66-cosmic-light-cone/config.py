"""Shared settings: units, the universes we draw, data locations and literature values.

Units throughout are **billions of years (Gyr)** for time and **billions of
light-years (Gly)** for distance. A light-year is the distance light covers in
one Julian year, so in these units the speed of light is exactly 1 Gly per Gyr:
a light ray on a spacetime diagram drawn in comoving distance against conformal
time is a 45-degree line.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
RAW_DIR = DATA_DIR / "raw"
LIT_DIR = DATA_DIR / "literature"
DERIVED_DIR = DATA_DIR / "derived"
FIG_DIR = HERE / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

ROOT = HERE.parent  # the atlas; sibling projects' caches are reused when present

# --- Units --------------------------------------------------------------------

C_KM_S = 299_792.458
#: 1 Mpc in Gly: 3.0856775814913673e22 m / (299792458 m/s x 365.25 x 86400 s x 1e9).
GLY_PER_MPC = 3.0856775814913673e22 / (299_792_458.0 * 365.25 * 86_400.0 * 1e9)
#: H0 in km/s/Mpc -> 1/Gyr.
KMSMPC_TO_PER_GYR = 1e3 / 3.0856775814913673e22 * 365.25 * 86_400.0 * 1e9


# --- The universes ------------------------------------------------------------


@dataclass(frozen=True)
class Params:
    """One homogeneous universe.

    ``Om`` is cold dark matter + baryons today, *excluding* massive neutrinos
    (astropy's ``Om0`` convention; Planck's quoted Omega_m includes them).
    Radiation is photons at ``T_cmb`` plus N_eff neutrino species; the species
    listed in ``m_nu`` have masses (eV) and each carries N_eff/3 of a species,
    the rest are massless (CAMB's and astropy's convention). ``T_cmb = 0`` turns
    radiation off, as in Davis & Lineweaver's (Omega_M, Omega_Lambda) models.
    Dark energy has w(a) = w0 + wa (1 - a) (a cosmological constant is -1, 0);
    its density today ``Ode`` closes the budget unless given, in which case
    curvature takes up the rest. ``w_future`` = "cpl" extrapolates w(a) to
    a > 1; "frozen" holds w at w0 after today.
    """

    name: str
    label: str
    H0: float                       # km/s/Mpc
    Om: float                       # CDM + baryons (no neutrinos)
    Ob: float = 0.0                 # baryons alone (only for CAMB comparisons)
    Ode: float | None = None        # None: by closure with Ok
    Ok: float = 0.0
    T_cmb: float = 2.7255
    N_eff: float = 3.046
    m_nu: tuple = (0.06,)
    w0: float = -1.0
    wa: float = 0.0
    w_future: str = "cpl"
    source: str = ""

    def but(self, **kw) -> "Params":
        return replace(self, **kw)


#: Planck 2018 TT,TE,EE+lowE+lensing+BAO (Planck Collaboration 2020, A&A 641, A6,
#: Table 2, last column) exactly as astropy's ``Planck18`` holds it: H0 = 67.66,
#: Om0 = 0.30966 (Omega_m 0.3111 with the 0.06 eV neutrino), Ob0 = 0.04897,
#: Tcmb 2.7255 K, N_eff 3.046.
PLANCK18 = Params("planck18", "Planck 2018", H0=67.66, Om=0.30966, Ob=0.04897,
                  source="Planck 2018 VI, Table 2 (TT,TE,EE+lowE+lensing+BAO) = astropy Planck18")

#: Davis & Lineweaver (2004), Fig. 1: (Omega_M, Omega_Lambda) = (0.3, 0.7), H0 = 70.
#: Whether they included radiation is tested against their figure (dl2004.py).
DL2004 = Params("dl2004", "Davis & Lineweaver 2004", H0=70.0, Om=0.3, Ode=0.7, T_cmb=0.0,
                N_eff=0.0, m_nu=(), source="Davis & Lineweaver 2004, PASA 21, 97, Fig. 1")
DL2004_RAD = DL2004.but(name="dl2004_rad", label="Davis & Lineweaver + radiation",
                        T_cmb=2.7255, N_eff=3.046, m_nu=(), Ode=None)

#: Planck's published values for the same column, used as checks.
PLANCK18_PUBLISHED = {
    "H0": (67.66, 0.42), "Om_incl_nu": (0.3111, 0.0056), "OL": (0.6889, 0.0056),
    "Om_h3": (0.09635, 0.00030), "age_gyr": (13.787, 0.020), "z_star": (1089.80, 0.21),
    "r_star_mpc": (144.57, 0.22), "z_drag": (1060.01, 0.29), "r_drag_mpc": (147.21, 0.23),
    "z_eq": (3387, 21),
}

#: The spec's numbers (root README, O2), all for Planck 2018.
SPEC_CHECKS = {"particle_horizon_gly": 46.2, "event_horizon_gly": 16.6, "z_hubble_sphere": 1.48,
               "z_min_angular_size": 1.59, "z_event_horizon": 1.85, "reachable_fraction": 0.05}

#: What Davis & Lineweaver (2004) state in their text for (0.3, 0.7), H0 = 70.
DL2004_TEXT = {
    "z_hubble_sphere": (1.46, "Sec. 2: 'all objects with redshift greater than z ~ 1.46 are receding faster than light'"),
    "z_event_horizon": (1.8, "Sec. 2: 'galaxies with redshift z ~ 1.8 are currently crossing our event horizon'"),
    "particle_horizon_gly": (46.0, "Fig. 1 caption: 'about 46 Glyr, the current distance to the particle horizon'"),
    "age_gyr": (13.5, "Fig. 1 caption: aspect ratio '46 Glyr/13.5 Gyr'"),
    "cmb_v_now_c": (3.2, "Sec. 3.3: CMB emitters' 'current recession velocity ... is 3.2c'"),
    "cmb_v_then_c": (58.1, "Sec. 3.3: 'At the time of emission their speed was 58.1c'"),
    "superluminal_emission_before_gyr": (5.0, "Sec. 3.3: photons 'emitted in the first ~five billion years' were emitted superluminally"),
    "symmetry_worldline_z": (10.0, "Sec. 3.4: particle and event horizons mirror each other about the z = 10 worldline"),
    "sr_sigma": (23.0, "Sec. 4.2: SR interpretation 'ruled out at more than 23 sigma'"),
    "linear_sigma": (12.0, "Fig. 5 caption: v = cz model 'ruled out at 12 sigma'"),
    "sr_dchi2": (530.0, "commented-out LaTeX in Fig. 5 caption: chi2 569 - 38.9"),
    "linear_dchi2": (154.0, "commented-out LaTeX in Fig. 5 caption: chi2 193 - 38.9"),
    "chi2_lcdm": (38.9, "commented-out LaTeX in Fig. 5 caption"),
    "script_M": (-3.45, "Sec. 4.2: 'We have taken M_B = -3.45'"),
}

#: Other universes for the sensitivity table, from this atlas's own fits.
OTHER_UNIVERSES = [
    # D3 (09-planck-cmb): our own CAMB fit to the Planck TT spectrum.
    Params("d3_tt", "D3: our Planck TT fit", H0=66.658, Om=0.32471 - 0.00141, Ob=0.0499873,
           source="09-planck-cmb/data/results.json fit (Omega_m 0.3247 incl. nu)"),
    # D5 (42-evolving-dark-energy): DESI DR2 BAO + compressed CMB, flat LCDM.
    Params("d5_lcdm", "D5: DESI DR2 + CMB, LCDM", H0=68.354, Om=0.30054 - 0.00138, Ob=0.0479680,
           source="42-evolving-dark-energy maps desi_cmb-desi lcdm"),
    # D1 (05-supernova-cosmology): SH0ES-calibrated H0 with the Pantheon+ flat Omega_m.
    Params("d1_shoes", "D1: SH0ES H0 + Pantheon+ flat", H0=73.55, Om=0.3315,
           source="05-supernova-cosmology: H0 73.55 +- 1.02, flat Om 0.332 +- 0.018"),
    # D5: DESI + CMB + DES-SN5YR w0wa best fit; w(a) extrapolated (CPL) or frozen at w0.
    Params("d5_w0wa_cpl", "D5: DESI+CMB+DESY5 w0wa, w(a) extrapolated", H0=66.732, Om=0.31868 - 0.00145,
           Ob=0.050058, w0=-0.76605, wa=-0.77818, w_future="cpl",
           source="42-evolving-dark-energy maps desi_cmb-desi_desy5 w0wa"),
    Params("d5_w0wa_frozen", "D5: same, w frozen at w0 after today", H0=66.732, Om=0.31868 - 0.00145,
           Ob=0.050058, w0=-0.76605, wa=-0.77818, w_future="frozen",
           source="42-evolving-dark-energy maps desi_cmb-desi_desy5 w0wa"),
]

#: Textbook universes for comparison curves.
EDS = Params("eds", "Einstein-de Sitter (matter only)", H0=67.66, Om=1.0, Ode=0.0, T_cmb=0.0,
             N_eff=0.0, m_nu=())
MILNE = Params("milne", "Milne (empty)", H0=67.66, Om=0.0, Ode=0.0, T_cmb=0.0, N_eff=0.0, m_nu=())
DE_SITTER = Params("de_sitter", "de Sitter (Lambda only)", H0=67.66, Om=0.0, Ode=1.0, T_cmb=0.0,
                   N_eff=0.0, m_nu=())

#: The worldlines Davis & Lineweaver label, by today's redshift of their light.
WORLDLINE_Z = (1.0, 3.0, 10.0, 1090.0)


# --- Data locations -------------------------------------------------------------

USER_AGENT = "astrophysics-atlas-O2/1.0 (educational; python-urllib)"

DL2004_ARXIV_SRC = "https://arxiv.org/e-print/astro-ph/0310808"   # v2, LaTeX + EPS figures
PLANCK18_ARXIV_PDF = "https://arxiv.org/pdf/1807.06209v4"

BAO_REPO_RAW = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/"
BAO_FILES = {
    "desi_dr2_mean.txt": "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
    "desi_dr2_cov.txt": "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
}
SN_REPO_RAW = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/"
SN_FILES = {"Pantheon+SH0ES.dat": "Pantheon%2BSH0ES.dat",
            "Pantheon+SH0ES_STAT+SYS.cov": "Pantheon%2BSH0ES_STAT%2BSYS.cov"}

GWOSC_URL = "https://gwosc.org/eventapi/json/GWTC/"

JADES_DR4_URL = "https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits"
JADES_README_URL = "https://jades.herts.ac.uk/DR4/Readme_DR4_catalogues.md"

VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
#: SDSS DR16Q (Lyke et al. 2020), redshift column only.
DR16Q_VIZIER = "VII/289/dr16q"
#: Farina et al. 2022, z > 5.7 quasars with Mg II and [C II] redshifts.
FARINA22_VIZIER = "J/ApJ/941/106/table1"

SKYSERVER_SQL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"

#: Sibling projects' caches, reused by hard link when present (same bytes, no download).
SIBLING_FILES = {
    "raw/bao/desi_dr2_mean.txt": "42-evolving-dark-energy/data/bao/desi_dr2_mean.txt",
    "raw/bao/desi_dr2_cov.txt": "42-evolving-dark-energy/data/bao/desi_dr2_cov.txt",
    "raw/sn/Pantheon+SH0ES.dat": "05-supernova-cosmology/data/Pantheon+SH0ES.dat",
    "raw/sn/Pantheon+SH0ES_STAT+SYS.cov": "05-supernova-cosmology/data/Pantheon+SH0ES_STAT+SYS.cov",
}
#: Larger derived tables from siblings, read in place when present (else fetched).
SIBLING_DR16Q = "34-quasars/data/derived/dr16q_columns.parquet"
SIBLING_MGS = "07-cosmic-web/data/mgs_dr17.npz"


# --- Landmark objects (literature redshifts, each with its reference) ----------


@dataclass(frozen=True)
class Landmark:
    key: str
    name: str
    kind: str
    z: float
    z_err: float
    ref: str
    note: str = ""
    label_side: str = "right"
    extra: dict = field(default_factory=dict)


#: Record-holders that are not in the catalogues we load. (JADES-GS-z14-0, GN-z11
#: and the other named JWST galaxies come from JADES DR4 itself, with its redshifts:
#: GS-z14-0 at 14.1796, the ALMA [O III] value of Carniani et al. 2025; Schouws et
#: al. 2025 find 14.1793 +- 0.0007 from the same line.)
LANDMARKS = [
    Landmark("mom_z14", "MoM-z14", "galaxy", 14.44, 0.02,
             "Naidu et al. 2026, Open J. Astrophys. 9 (arXiv:2505.11263)",
             "most distant spectroscopically confirmed galaxy (rest-UV lines, NIRSpec prism)"),
    Landmark("j0313", "J0313-1806", "quasar", 7.6423, 0.0013, "Wang et al. 2021, ApJL 907, L1",
             "most distant known quasar ([C II] redshift)"),
    Landmark("grb090423", "GRB 090423", "grb", 8.23, 0.07,
             "Tanvir et al. 2009, Nature 461, 1254; Salvaterra et al. 2009, Nature 461, 1258",
             "most distant spectroscopic gamma-ray burst"),
    Landmark("grb250314a", "GRB 250314A + SN", "grb", 7.3, 0.1,
             "Levan et al. 2025, A&A 704, L8; Cordier et al. 2025, A&A 704, L7",
             "most distant supernova (seen by JWST after the burst)"),
    Landmark("sn2023adsy", "SN 2023adsy", "supernova", 2.903, 0.007, "Pierel et al. 2024, ApJL 971, L32",
             "most distant spectroscopic Type Ia supernova"),
]
