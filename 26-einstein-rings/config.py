"""Project configuration: paths, data sources, conventions and every assumed number.

Everything the analysis *assumes* rather than measures lives in this file, so
the assumptions can be read, and changed, in one place. Most of them are the
conventions of the SLACS papers we compare against, kept so that a difference
between our numbers and theirs means a difference in the measurement, not in
the bookkeeping.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
TABLE_DIR = DATA_DIR / "tables"  # SLACS catalogues from CDS
LIT_DIR = DATA_DIR / "literature"  # published per-lens results, for grading only
SDSS_DIR = DATA_DIR / "sdss"  # SDSS DR17 spectra of every SLACS candidate
HST_DIR = DATA_DIR / "hst"  # Hubble cutouts of every grade-A lens
MODEL_DIR = DATA_DIR / "lens_models"  # our own lens-model fits, cached
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, TABLE_DIR, LIT_DIR, SDSS_DIR, HST_DIR, MODEL_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-project-atlas/J5 (hobby research; python-urllib)"


# --------------------------------------------------------------------------
# Data sources
# --------------------------------------------------------------------------

CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp/"
BOLTON08 = "Bolton et al. 2008, ApJ 682, 964 (SLACS V)"
AUGER09 = "Auger et al. 2009, ApJ 705, 1099 (SLACS IX)"
AUGER10 = "Auger et al. 2010, ApJ 724, 511 (SLACS X)"
TREU10 = "Treu et al. 2010, ApJ 709, 1195"


@dataclass(frozen=True)
class Source:
    """One file we download, and why."""

    #: Path under ``data/``.
    filename: str
    url: str
    #: What it holds and what we use it for.
    what: str
    #: The paper or archive to cite for it.
    reference: str
    #: ``True`` if the analysis is built on it; ``False`` if it is a published
    #: result used only to grade our own at the end.
    input: bool = True
    #: Payload check: "cds", "fits", "gzip" or "zip".
    kind: str = "cds"


def _cds(catalog: str, name: str, what: str, reference: str) -> Source:
    folder = catalog.replace("/", "_")
    return Source(f"tables/{folder}/{name}", f"{CDS_FTP}{catalog}/{name}", what, reference)


TABLE_SOURCES: list[Source] = [
    _cds("J/ApJ/682/964", "ReadMe", "Column definitions for the SLACS V tables", BOLTON08),
    _cds("J/ApJ/682/964", "table4.dat",
         "All 131 SLACS candidates observed with HST: SDSS plate-MJD-fibre, both redshifts, "
         "F814W de Vaucouleurs magnitude and effective radius, velocity dispersion, "
         "morphology and lens grade (A/B/X)", BOLTON08),
    _cds("J/ApJ/682/964", "table5.dat",
         "Singular-isothermal-ellipsoid lens models of 63 grade-A lenses: Einstein radius, "
         "axis ratio, position angle; ring-subset and good-sigma flags", BOLTON08),
    _cds("J/ApJ/682/964", "notes.dat", "Why some grade-A lenses were not modelled", BOLTON08),
    _cds("J/ApJ/705/1099", "ReadMe", "Column definitions for the SLACS IX tables", AUGER09),
    _cds("J/ApJ/705/1099", "table1.dat", "12 SLACS targets added after SLACS V", AUGER09),
    _cds("J/ApJ/705/1099", "table3.dat",
         "85 grade-A lenses: redshifts, sigma, morphology, HST magnitudes and effective "
         "radii in B, V, I, H", AUGER09),
    _cds("J/ApJ/705/1099", "table4.dat",
         "85 grade-A lenses: Einstein radius in kpc, mass inside it, stellar masses for "
         "Chabrier and Salpeter IMFs, and their stellar fraction inside the ring", AUGER09),
]

#: arXiv LaTeX sources of two papers whose per-lens tables are not on CDS. The
#: analysis never reads them; they grade our density slopes and IMF numbers.
LITERATURE_SOURCES: list[Source] = [
    Source("literature/auger2010_1007.2880.tar.gz", "https://arxiv.org/e-print/1007.2880",
           "LaTeX source; Table 1 lists each lens's density slope gamma' and dark fraction "
           "inside r_e/2", AUGER10, input=False, kind="gzip"),
    Source("literature/treu2010_0911.3392.tar.gz", "https://arxiv.org/e-print/0911.3392",
           "LaTeX source; Table 1 lists each lens's IMF mismatch parameter alpha", TREU10,
           input=False, kind="gzip"),
]

#: SDSS DR17 spectra ("lite" format: the coadded spectrum and the pipeline's
#: best-fit galaxy model). Every SLACS candidate is a DR7 legacy plate, run 26.
SDSS_SPECTRUM_URL = ("https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/"
                     "{plate:04d}/spec-{plate:04d}-{mjd:05d}-{fiber:04d}.fits")
SDSS_REFERENCE = "SDSS DR17 (Abdurro'uf et al. 2022); SDSS-I/II legacy spectra (Abazajian et al. 2009)"

#: MAST's HST cutout service. It cuts every Hubble Advanced Products
#: multi-visit mosaic (MVM) that covers the position -- one FITS file per
#: instrument and filter, each with SCI (e-/s) and WHT (exposure weight)
#: extensions, drizzled to 0.04"/pixel -- and returns them as a zip.
HAPCUT_URL = "https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={ra:.6f}&dec={dec:.6f}&x={n}&y={n}&units=px"
HAPCUT_REFERENCE = ("HST programs 10174, 10494, 10587, 10798, 10886 (SLACS; PIs Koopmans, "
                    "Bolton) and others, via MAST Hubble Advanced Products")
HST_CUTOUT_PIXELS = 200  # 8" on a side at 0.04"/pixel
HST_PIXEL_ARCSEC = 0.04


# --------------------------------------------------------------------------
# Cosmology
# --------------------------------------------------------------------------

#: The SLACS papers' cosmology. We use it whenever we compare with them, so
#: that kpc and solar masses mean the same thing on both sides.
SLACS_H0 = 70.0
SLACS_OMEGA_M = 0.3

#: Planck 2018 (TT,TE,EE+lowE+lensing), for showing how little the choice
#: matters. Project D3 measured H0 = 66.7 +/- 0.9 from the same satellite.
PLANCK_H0 = 67.36
PLANCK_OMEGA_M = 0.3153


# --------------------------------------------------------------------------
# Velocity dispersions and dynamics (Koopmans et al. 2006, 2009 conventions)
# --------------------------------------------------------------------------

#: The SDSS spectrograph fibre has a 3" diameter.
FIBRE_RADIUS_ARCSEC = 1.5
#: Seeing during the SDSS spectroscopy, as assumed by Koopmans et al. 2006/2009.
#: Bolton et al. (2008) quote a median of ~1.4". With a 3" fibre it barely matters.
SEEING_FWHM_ARCSEC = 1.5
#: Aperture correction for velocity dispersions (Jorgensen et al. 1995):
#: sigma_ap / sigma_e8 = (R_ap / (R_e/8))^-0.04.
JORGENSEN_SLOPE = -0.04
#: Bolton et al. (2008b) weighted their sigma averages with errors floored at 5%.
SIGMA_ERROR_FLOOR = 0.05
#: A Hernquist profile with a = R_e / 1.8153 has the projected half-light
#: radius R_e (Koopmans et al. 2006).
HERNQUIST_RE_OVER_A = 1.8153
#: Range of total-density slopes searched (Koopmans et al. 2009 used 1.1-2.9).
GAMMA_RANGE = (1.1, 2.9)
#: Treu et al. (2010): NFW dark haloes with the break radius fixed at 30 kpc,
#: the expected value for the ~1e13 solar-mass haloes of SLACS lenses.
NFW_BREAK_RADIUS_KPC = 30.0

#: Systematic floor on Einstein radii adopted by Bolton et al. (2008): 2%.
EINSTEIN_RADIUS_ERROR = 0.02


# --------------------------------------------------------------------------
# Our own lens models
# --------------------------------------------------------------------------

#: Model the central part of each cutout: 6.4" square.
MODEL_BOX_PIXELS = 160
#: F814W point-spread function, approximated as a circular Moffat profile.
#: Drizzled ACS/WFC F814W stars have FWHM ~0.10" (ACS instrument handbook).
PSF_FWHM_ARCSEC = 0.10
PSF_MOFFAT_BETA = 2.5
#: Rays per pixel side when evaluating the lensed source.
SUPERSAMPLE = 2

#: Lenses shown individually, chosen for variety before any comparison with the
#: published models; two whose fits flagged themselves (source size at its
#: bound) were then swapped for lenses with unflagged fits.
SHOWCASE = [
    "J0946+1006",  # the "Jackpot": two sources at two redshifts, two rings
    "J1627-0053",  # a near-complete ring
    "J0737+3216",
    "J1430+4105",
    "J2300+0022",
    "J0912+0029",
    "J1630+4520",  # the largest Einstein radius in SLACS
    "J2321-0939",  # one of the nearest lenses, z = 0.08
    "J0029-0055",
    "J0252+0039",
    "J1205+4910",  # one where our model and SLACS's disagree (see section 4)
    "J0037-0942",  # a pair of images rather than a ring
]

#: The lens whose SDSS spectrum demonstrates two redshifts in one fibre.
SPECTRUM_LENS = "J1205+4910"

#: Neighbouring galaxies and stars masked by hand before our lens fits, as
#: (arcsec east, arcsec north, radius in arcsec) from the lens centre. Chosen by
#: eye from the images; see the README for when and why each was added.
#
#: The rule, applied to every lens on a blind inspection of the lens-subtracted
#: images (our own fits only, no published values shown): mask objects that are
#: clearly separate -- compact sources, or elongated ones not curved around the
#: lens -- and never anything arc-like centred on it. The headline comparison
#: uses the fits *without* these masks; the masked fits are reported alongside.
MANUAL_MASKS: dict[str, list[tuple[float, float, float]]] = {
    # a broad band ~1.2" south of the lens, 5" long, curved about a point far from
    # the lens: a neighbouring galaxy or tidal debris, not a ring around this lens
    "J0157-0056": [(2.6, -1.25, 0.35), (2.0, -1.2, 0.35), (1.4, -1.2, 0.35), (0.8, -1.25, 0.3),
                   (0.2, -1.2, 0.3), (-0.4, -1.05, 0.3), (-1.0, -0.9, 0.3), (-1.6, -0.75, 0.3)],
    # compact source 2.35" west
    "J0728+3835": [(-2.35, 0.1, 0.3)],
    # neighbour galaxy 2" north-west, and a compact blob 1.8" north-east
    "J0808+4706": [(-0.75, 2.1, 0.4), (-0.65, 1.65, 0.4), (1.0, 1.5, 0.25)],
    # a second galaxy, about as bright as the lens, 1" south
    "J1618+4353": [(0.0, -1.0, 0.45)],
    # a resolved neighbour galaxy 1.6" south-west
    "J1718+6424": [(-1.25, -1.05, 0.45)],
}
