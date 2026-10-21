"""Project configuration: paths, the data sources, and the published numbers we
grade ourselves against.

Two independent teams have followed the star S2 around Sagittarius A* since the
1990s, with different telescopes, cameras and reference frames:

* **UCLA / Keck** (Andrea Ghez's group), who call the star S0-2. Do et al.
  (2019, Science 365, 664) released every position and radial velocity behind
  their gravitational-redshift paper, 1995-2019, together with their posterior
  chains. The group's old web page is gone; the file survives in the Internet
  Archive, which is where we fetch it from (the SHA-256 in the manifest pins
  the exact bytes).
* **MPE / ESO VLT** (Reinhard Genzel's group), who call it S2. Gillessen et al.
  (2017, ApJ 837, 30) published their 1992-2016 data for 17 stars, S2
  included, through the CDS VizieR service.

The two data sets share nothing but the sky, so each can weigh the black hole
on its own; fitting them together needs a separate reference frame for each.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CHAIN_DIR = DATA_DIR / "chains"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"
VIEWER_PATH = FIGURE_DIR / "s2_orbit_explorer.html"

for _d in (DATA_DIR, RAW_DIR, CHAIN_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Data sources
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Source:
    """One file we download, and why."""

    #: File name under ``data/raw/``.
    filename: str
    url: str
    #: What it holds and what we use it for.
    what: str
    #: The paper to cite for it.
    reference: str
    #: Where the file originally lived, if we fetch an archived copy.
    original_url: str = ""
    #: A string the payload must contain (guards against HTML error pages).
    marker: bytes = b""


DO2019_PAPER = ("Do, T., Hees, A., Ghez, A. et al. 2019, Science 365, 664 "
                "(arXiv:1907.10731), 'Relativistic redshift of the star S0-2 orbiting "
                "the Galactic center supermassive black hole'")
GILLESSEN2017_PAPER = ("Gillessen, S., Plewa, P. M., Eisenhauer, F. et al. 2017, ApJ 837, 30 "
                       "(arXiv:1611.09144), 'An update on monitoring stellar orbits in the "
                       "Galactic Center'; VizieR J/ApJ/837/30, doi:10.26093/cds/vizier.18370030")

DO2019_ORIGINAL = "http://www.astro.ucla.edu/~ghezgroup/gc/gr/data.zip"
VIZIER_G17 = "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/837/30/"

SOURCES: list[Source] = [
    Source(
        "do2019_data.zip",
        "https://web.archive.org/web/20190901075942id_/" + DO2019_ORIGINAL,
        "UCLA/Keck release for Do et al. 2019: S0-2_astro.csv (46 Keck positions, "
        "1995-2018), S0-2_rv.csv (116 radial velocities from Keck, Subaru, Gemini and "
        "the VLT, 2000-2019) and posterior_chains_data.zip (the authors' MultiNest "
        "chains for their Newtonian, GR and free-redshift fits)",
        DO2019_PAPER,
        original_url=DO2019_ORIGINAL,
        marker=b"S0-2_astro.csv",
    ),
    Source(
        "gillessen2017_ReadMe.txt", VIZIER_G17 + "ReadMe",
        "Column definitions for the two tables below",
        GILLESSEN2017_PAPER, marker=b"J/ApJ/837/30",
    ),
    Source(
        "gillessen2017_table5.dat", VIZIER_G17 + "table5.dat",
        "MPE astrometry (NTT 1992-2002, VLT/NACO 2002-2016) and radial velocities "
        "(VLT/SINFONI, plus three Keck points) for 17 S-stars; we use S2's 145 "
        "positions and 44 velocities",
        GILLESSEN2017_PAPER, marker=b"1992.224",
    ),
    Source(
        "gillessen2017_table3.dat", VIZIER_G17 + "table3.dat",
        "Published orbital elements of the 40 S-stars with determined orbits: used for "
        "the 'every star feels the same mass' figure, not for the S2 fit",
        GILLESSEN2017_PAPER, marker=b"S2 ",
    ),
]


# --------------------------------------------------------------------------
# Data selection
# --------------------------------------------------------------------------

#: Do et al. (2019, supplement section 1.2.5) exclude this Keck epoch because
#: Sgr A* was unusually bright that night and likely biased S0-2's position.
#: The released CSV still contains it; we drop it the same way.
KECK_EXCLUDED_EPOCHS = [2018.673]

#: Gillessen et al. 2017 (table 1, footnote): the NACO astrometry is already
#: down-weighted by a global factor 1.42, and S2's 2002 positions (confused
#: near pericentre) carry inflated errors. We use the released errors as they
#: are and let an overall scale factor absorb any remaining mismatch.


# --------------------------------------------------------------------------
# Published numbers we grade ourselves against (never used in a fit)
# --------------------------------------------------------------------------

PUBLISHED = {
    # Do et al. 2019, Table S3 ("estimation" column) and main text.
    "do2019_redshift": {"gm6": (3.975, 0.058, 0.026), "r0": (7.959, 0.059, 0.032),
                        "f": (0.88, 0.16, 0.047),
                        "ref": "Do et al. 2019, Science 365, 664 (Y free)"},
    "do2019_gr": {"gm6": (3.964, 0.047, 0.026), "r0": (7.946, 0.050, 0.032),
                  "ref": "Do et al. 2019 (GR assumed)"},
    # Gillessen et al. 2017, Table 1.
    "g17_s2_vlt_nopriors": {"r0": (8.17, 0.20), "gm6": (4.25, 0.20),
                            "alpha_mas": (0.23, 0.39), "delta_mas": (-2.10, 0.61),
                            "valpha_uas_yr": (88, 40), "vdelta_uas_yr": (-2, 63),
                            "vz_kms": (28.3, 7.0), "reduced_chi2": 1.19,
                            "ref": "Gillessen et al. 2017, Table 1 row 1 (S2, VLT data, no priors)"},
    "g17_s2_combined_nopriors": {"r0": (8.17, 0.15), "gm6": (4.30, 0.15),
                                 "ref": "Gillessen et al. 2017, Table 1 row 4 (S2, VLT + Keck, no priors)"},
    "g17_multistar": {"r0": (8.32, 0.07, 0.14), "gm6": (4.28, 0.10, 0.21),
                      "ref": "Gillessen et al. 2017 (17 stars, with priors)"},
    "g17_s2_elements": {"a_arcsec": (0.1255, 0.0009), "e": (0.8839, 0.0019),
                        "inc": (134.18, 0.40), "Omega": (226.94, 0.60), "omega": (65.51, 0.57),
                        "tp": (2002.33, 0.01), "period": (16.00, 0.02),
                        "ref": "Gillessen et al. 2017, Table 3"},
    # GRAVITY Collaboration papers (VLTI interferometry, 2017 onwards: data not public).
    "gravity2018": {"gm6": (4.100, 0.034), "r0": (8.122, 0.031), "f": (0.90, 0.09, 0.15),
                    "ref": "GRAVITY Collaboration 2018, A&A 615, L15"},
    "gravity2019": {"gm6": (4.152, 0.014), "r0": (8.178, 0.013, 0.022),
                    "ref": "GRAVITY Collaboration 2019, A&A 625, L10"},
    "gravity2020": {"gm6": (4.261, 0.012), "r0": (8.246, 0.009), "f_sp": (1.10, 0.19),
                    "ref": "GRAVITY Collaboration 2020, A&A 636, L5"},
    "gravity2022": {"gm6": (4.297, 0.012), "r0": (8.277, 0.009),
                    "ref": "GRAVITY Collaboration 2022, A&A 657, L12"},
    # The Event Horizon Telescope's image of Sgr A* (for the shadow prediction).
    "eht2022": {"ring_uas": (51.8, 2.3), "shadow_uas": (48.7, 7.0),
                "theta_g_uas": (4.8, 0.7, 1.4),
                "ref": "Event Horizon Telescope Collaboration 2022, ApJL 930, L12"},
}

#: The 2020 Nobel Prize in Physics: half to Roger Penrose, half jointly to
#: Reinhard Genzel and Andrea Ghez "for the discovery of a supermassive compact
#: object at the centre of our galaxy".
NOBEL_2020 = "Nobel Prize in Physics 2020 (Genzel and Ghez, with Penrose)"
