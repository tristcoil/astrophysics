"""Project configuration: paths, data sources, the galaxies, and every assumption.

Every number the project *chooses* rather than measures lives in this file: the
display stretch, which filter goes to which colour channel, what counts as
"white", the vision thresholds. A colour picture of a galaxy is a stack of such
choices, and the point of this project is to make each one visible and
changeable in one place.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlencode

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
CUTOUT_DIR = DATA_DIR / "cutouts"
CATALOG_DIR = DATA_DIR / "catalogs"
REFERENCE_DIR = DATA_DIR / "reference"
MULTI_DIR = DATA_DIR / "multiwavelength"
FIGURE_DIR = PROJECT_DIR / "figures"
GALLERY_DIR = FIGURE_DIR / "gallery"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, CUTOUT_DIR, CATALOG_DIR, REFERENCE_DIR, MULTI_DIR, FIGURE_DIR, GALLERY_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# The imaging survey
# --------------------------------------------------------------------------

#: The Legacy Surveys sky viewer (Dey et al. 2019, AJ 157, 168) serves FITS and
#: JPEG cutouts of several surveys on any grid we ask for. We use its "sdss"
#: layer: Sloan Digital Sky Survey DR8 imaging frames, coadded by the Legacy
#: Surveys team onto a common pixel grid ("sdssco" bricks, 0.396"/pixel) and
#: resampled to the requested grid on the fly.
VIEWER_URL = "https://www.legacysurvey.org/viewer"
SURVEY_LAYER = "sdss"

#: SDSS native pixel scale, arcsec. The cutout service resamples but keeps the
#: pixel values in nanomaggies *per native pixel* -- i.e. it preserves surface
#: brightness, not flux. Verified in ``analysis.validate_star_photometry``
#: against catalogue photometry of hundreds of stars.
NATIVE_PIXSCALE = 0.396

#: The three filters the SDSS layer provides, blue to red.
BANDS = ("g", "r", "i")

#: SDSS magnitudes are *nearly* on the AB system. The survey documents offsets
#: of u_AB = u - 0.04 and z_AB = z + 0.02, and says g, r, i are AB to within
#: about 0.01 mag -- the level of the whole-survey calibration itself. We treat
#: g, r, i as AB.
SDSS_AB_NOTE = "g, r, i treated as AB (documented offsets < 0.01 mag)"

#: The nanomaggy: flux such that m_AB = 22.5 - 2.5 log10(flux). One nanomaggy
#: is 3.631 microjansky.
NANOMAGGY_JY = 3.631e-6


# --------------------------------------------------------------------------
# The gallery
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Galaxy:
    """One galaxy we fetch, composite and analyse."""

    #: Short key used in file names, e.g. "m51".
    key: str
    #: Display name, e.g. "M51".
    name: str
    #: Popular name.
    common: str
    #: Catalogue designation used to look it up in RC3, e.g. "NGC 5194".
    ngc: str
    #: Centre of the cutout, degrees (J2000). Galaxy positions from SIMBAD.
    ra: float
    dec: float
    #: Side of the square cutout, arcmin.
    fov_arcmin: float
    #: Output pixel scale, arcsec/pixel.
    pixscale: float
    #: What kind of galaxy it is, in a few words.
    kind: str
    #: What it shows about colour.
    why: str
    #: Names of its integrated spectra in the Brown et al. (2014) atlas, if any.
    brown: tuple[str, ...] = ()
    #: Other catalogue galaxies inside the frame (RC3 names), for masks and labels.
    companions: tuple[str, ...] = ()

    @property
    def size_px(self) -> int:
        return int(round(self.fov_arcmin * 60.0 / self.pixscale))


#: Twelve galaxies spanning what galaxies look like, in the gallery's order:
#: face-on spirals; bulges and dust; edge-on discs; a starburst; a merger; and
#: a giant elliptical. All lie inside the SDSS footprint, so every panel
#: of the gallery uses the same camera, the same three filters and the same
#: calibration -- colours mean the same thing in every panel. Nine have an
#: integrated spectrum in the Brown et al. (2014) atlas, which gives an
#: independent check on the colours we measure from the images.
GALAXIES: list[Galaxy] = [
    Galaxy("m51", "M51", "Whirlpool Galaxy", "NGC 5194", 202.4840, 47.2050, 16.0, 0.70,
           "grand-design spiral and its companion",
           "Blue spiral arms of young stars, brown dust lanes, and a yellow companion "
           "(NGC 5195) made of old stars, caught mid-encounter.",
           brown=("NGC_5194", "NGC_5195"), companions=("NGC 5195",)),
    Galaxy("m101", "M101", "Pinwheel Galaxy", "NGC 5457", 210.80243, 54.34875, 34.0, 1.50,
           "giant face-on spiral",
           "A lopsided disc studded with giant star-forming regions; the faint outer arms "
           "test how deep the stretch reaches."),
    Galaxy("m74", "M74", "Phantom Galaxy", "NGC 628", 24.17394, 15.78364, 13.0, 0.60,
           "face-on grand-design spiral",
           "The textbook spiral: seen face-on, with a yellowish centre and blue arms.",
           brown=("NGC_0628",)),
    Galaxy("m63", "M63", "Sunflower Galaxy", "NGC 5055", 198.95531, 42.02937, 15.0, 0.70,
           "flocculent spiral",
           "No grand spiral pattern: short, patchy arm segments of young stars.",
           brown=("NGC_5055",)),
    Galaxy("m81", "M81", "Bode's Galaxy", "NGC 3031", 148.88822, 69.06530, 30.0, 1.30,
           "spiral with a large bulge",
           "A big, smooth, yellow bulge of old stars wrapped in thin blue arms."),
    Galaxy("m64", "M64", "Black Eye Galaxy", "NGC 4826", 194.18207, 21.68266, 12.0, 0.55,
           "spiral with a dark dust band",
           "A band of dust across the bright inner disc blocks and reddens the starlight "
           "behind it.",
           brown=("NGC_4826",)),
    Galaxy("m104", "M104", "Sombrero Galaxy", "NGC 4594", 189.99763, -11.62305, 11.0, 0.50,
           "edge-on spiral with a huge bulge",
           "A thin dust lane seen edge-on against an enormous glow of old stars.",
           brown=("NGC_4594",)),
    Galaxy("ngc5866", "NGC 5866", "Spindle Galaxy", "NGC 5866", 226.62317, 55.76331, 7.0, 0.40,
           "edge-on lenticular",
           "A disc galaxy that has stopped forming stars: old, yellow, with one razor-thin "
           "dust lane and no spiral arms.",
           brown=("NGC_5866",)),
    Galaxy("ngc4631", "NGC 4631", "Whale Galaxy", "NGC 4631", 190.53337, 32.54151, 17.0, 0.75,
           "edge-on star-forming disc",
           "A disc seen exactly edge-on, forming stars along its whole length; its small "
           "companion NGC 4627 sits just above it.",
           brown=("NGC_4631",), companions=("NGC 4627",)),
    Galaxy("m82", "M82", "Cigar Galaxy", "NGC 3034", 148.96846, 69.67970, 13.0, 0.60,
           "starburst galaxy",
           "Forming stars ten times faster than the Milky Way; its disc is choked with "
           "dust, and hydrogen gas blown out of it glows in the r filter."),
    Galaxy("mice", "The Mice", "NGC 4676", "NGC 4676A", 191.54440, 30.72730, 5.0, 0.40,
           "colliding pair",
           "Two spirals 300 million light-years away in the middle of merging, trailing "
           "long tidal tails of stars torn out by the collision. Measurements refer to the "
           "northern member, NGC 4676A, with its partner B treated as a companion.",
           brown=("NGC_4676_A",), companions=("NGC 4676B",)),
    Galaxy("m87", "M87", "Virgo A", "NGC 4486", 187.70593, 12.39112, 12.0, 0.55,
           "giant elliptical",
           "Almost only old stars: one smooth yellow-white glow, surrounded by thousands "
           "of globular clusters that show up as faint dots.",
           brown=("NGC_4486",)),
]

#: The galaxy the explanatory figures are built around.
FLAGSHIP = "m51"


def galaxy_by_key(key: str) -> Galaxy:
    for g in GALAXIES:
        if g.key == key:
            return g
    raise KeyError(f"{key!r} is not a gallery galaxy; known: {[g.key for g in GALAXIES]}")


# --------------------------------------------------------------------------
# Composite choices (all of them)
# --------------------------------------------------------------------------

#: Which filter drives which display channel. "Chromatic order": the shortest
#: wavelength goes to blue, the longest to red, as in every SDSS colour image.
#: It is *not* what the eye would see -- i (763 nm) is almost invisible to the
#: eye -- which is what ``color_physics`` and figure 06 are about.
CHANNEL_OF_BAND = {"i": "R", "r": "G", "g": "B"}

#: What "white" means. Each band is multiplied by a weight such that a star
#: with the Sun's colours comes out exactly neutral. The weights are computed
#: from the measured solar spectrum (CALSPEC) through the SDSS filter curves, so
#: this is a physical statement: grey in these pictures means "the colour of
#: sunlight"; bluer than the Sun shows blue, redder shows yellow to red.
WHITE_POINT = "sun"

#: The stretch, identical for every galaxy and expressed in physical units, so
#: that equal surface brightness gives equal pixel brightness in every panel.
#: This is the arcsinh stretch of Lupton et al. (2004): linear for faint light,
#: logarithmic for bright light, applied to the *total* intensity so that each
#: pixel keeps its colour. The stretch below maps the sky (mu_r ~ 25
#: mag/arcsec^2) to near black and mu_r ~ 16.5 (the brightest galaxy cores) to
#: full brightness. The intensity unit is nanomaggies per arcsec^2 in r.
LUPTON_Q = 8.0
#: Surface brightness (nanomaggies/arcsec^2, solar-weighted mean of the three
#: channels) at which the display saturates.
LUPTON_SATURATION = 250.0

#: The Legacy Surveys viewer's own JPEG algorithm for the SDSS layer (legacypipe
#: ``sdss_rgb`` defaults): per-band scales, an offset m, and Q. We re-implement
#: it only to prove we read the data exactly as the survey does.
VIEWER_SDSS_RGB = {"scales": {"g": (2, 2.5), "r": (1, 1.5), "i": (0, 1.0)}, "m": 0.02, "Q": 20.0}


# --------------------------------------------------------------------------
# Vision
# --------------------------------------------------------------------------

#: Maximum luminous efficacy, lm/W: photopic (cone) and scotopic (rod) vision.
K_M_PHOTOPIC = 683.0
K_M_SCOTOPIC = 1700.0

#: CIE 191:2010 "Recommended system for mesopic photometry": vision is
#: *mesopic* -- rods and cones both active, colour vision degraded -- between
#: these photopic luminances. Below the lower limit only rods work and there is
#: no colour at all; above the upper one vision is fully photopic.
MESOPIC_LOWER_CD_M2 = 0.005
MESOPIC_UPPER_CD_M2 = 5.0

#: Everyday scenes for context, each defined by the numbers that produce its
#: luminance so nothing is quoted from memory: a white sheet of paper (a matte,
#: "Lambertian" surface reflecting 80% of the light) has luminance
#: 0.8 x illuminance / pi. The natural moonless sky is 22.0 V mag/arcsec^2, the
#: value Falchi et al. (2016, Sci. Adv. 2, e1600377) adopt for an unpolluted
#: sky; it is converted with the same formula the galaxies go through.
PAPER_REFLECTANCE = 0.8
CONTEXT_SCENES = {
    "white paper in full sunlight": {"illuminance_lux": 100_000.0},
    "white paper on an office desk": {"illuminance_lux": 500.0},
    "white paper by the light of a full Moon": {"illuminance_lux": 0.25},
    "the natural night sky, no Moon": {"v_mag_arcsec2": 22.0},
}


# --------------------------------------------------------------------------
# Data sources
# --------------------------------------------------------------------------

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
    #: What kind of payload to expect ("fits", "jpeg", "csv", "text", "tsv", "xml").
    kind: str = "fits"
    #: Seconds to wait before this request (politeness towards small services).
    pause: float = 0.0


def cutout_url(g: Galaxy, *, layer: str = SURVEY_LAYER, kind: str = "fits",
               pixscale: float | None = None, size: int | None = None) -> str:
    pixscale = pixscale or g.pixscale
    size = size or g.size_px
    endpoint = "fits-cutout" if kind == "fits" else "jpeg-cutout"
    return f"{VIEWER_URL}/{endpoint}?" + urlencode(
        {"ra": f"{g.ra:.5f}", "dec": f"{g.dec:.5f}", "layer": layer,
         "pixscale": f"{pixscale:.3f}", "size": size})


#: The viewer cuts off any request that takes longer than about two minutes,
#: and large SDSS cutouts can take longer than that to render and send. So each
#: galaxy is fetched as a grid of tiles no bigger than this (about 2.5 MB each),
#: every tile a little larger than its share of the frame, and the tiles are
#: stitched onto the exact target grid with ``reproject`` (``imaging.py``).
TILE_MAX_PX = 460
TILE_MARGIN_PX = 12


@dataclass(frozen=True)
class Tile:
    """One piece of a galaxy's frame: the region [x0, x1) x [y0, y1) of the
    target grid it supplies, and the (larger) cutout requested to cover it."""

    key: str
    ix: int
    iy: int
    x0: int
    x1: int
    y0: int
    y1: int
    ra: float
    dec: float
    size: int

    @property
    def filename(self) -> str:
        return f"cutouts/{self.key}/tile_{self.iy}_{self.ix}.fits"


def target_wcs(g: Galaxy):
    """The frame every image of this galaxy is put on: a gnomonic (TAN)
    projection centred on the galaxy, north up and east left -- the same
    convention the viewer uses for its cutouts."""
    from astropy.wcs import WCS

    n = g.size_px
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [g.ra, g.dec]
    w.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
    w.wcs.cd = [[-g.pixscale / 3600.0, 0.0], [0.0, g.pixscale / 3600.0]]
    return w


def tile_layout(g: Galaxy) -> list[Tile]:
    import numpy as np

    n = g.size_px
    k = int(np.ceil(n / TILE_MAX_PX))
    edges = np.linspace(0, n, k + 1).round().astype(int)
    core = int(np.max(np.diff(edges)))
    size = core + 2 * TILE_MARGIN_PX
    w = target_wcs(g)
    tiles = []
    for iy in range(k):
        for ix in range(k):
            x0, x1, y0, y1 = edges[ix], edges[ix + 1], edges[iy], edges[iy + 1]
            cx, cy = (x0 + x1 - 1) / 2.0, (y0 + y1 - 1) / 2.0
            ra, dec = w.pixel_to_world_values(cx, cy)
            tiles.append(Tile(g.key, ix, iy, int(x0), int(x1), int(y0), int(y1),
                              float(ra), float(dec), size))
    return tiles


def central_tile(g: Galaxy) -> Tile:
    tiles = tile_layout(g)
    k = int(round(len(tiles) ** 0.5))
    return next(t for t in tiles if t.ix == k // 2 and t.iy == k // 2)


def tile_url(t: Tile, pixscale: float, *, kind: str = "fits") -> str:
    endpoint = "fits-cutout" if kind == "fits" else "jpeg-cutout"
    return f"{VIEWER_URL}/{endpoint}?" + urlencode(
        {"ra": f"{t.ra:.6f}", "dec": f"{t.dec:.6f}", "layer": SURVEY_LAYER,
         "pixscale": f"{pixscale:.3f}", "size": t.size})


#: SDSS DR17 SkyServer SQL: clean point sources with PSF magnitudes, for the
#: photometric check of the cutouts.
SKYSERVER_SQL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"

#: Gaia DR3 TAP. Gaia sources for masking foreground stars, joined with the
#: Gaia Synthetic Photometry Catalogue (GSPC; Gaia Collaboration, Montegriffo et
#: al. 2023), which gives SDSS-system magnitudes computed from Gaia's own
#: spectrophotometry -- a calibration independent of SDSS.
GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"

VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"

#: IRSA's Galactic dust reddening service: E(B-V) through the Milky Way towards
#: each galaxy, from the Schlegel, Finkbeiner & Davis (1998) dust map.
IRSA_DUST = "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust"

#: Extinction per unit E(B-V)_SFD in the SDSS bands for a Milky-Way dust law
#: (R_V = 3.1): Schlafly & Finkbeiner 2011, ApJ 737, 103, Table 6. These already
#: include their 14% recalibration of the SFD map, so they multiply the
#: *original* SFD value.
EXTINCTION_PER_EBV_SFD = {"u": 4.239, "g": 3.303, "r": 2.285, "i": 1.698, "z": 1.263}


def _sdss_star_query(g: Galaxy) -> str:
    half = g.fov_arcmin / 2 / 60.0
    import math
    dra = half / max(math.cos(math.radians(g.dec)), 0.1)
    sql = (
        "SELECT p.objid, p.ra, p.dec, p.psfMag_g, p.psfMag_r, p.psfMag_i, "
        "p.psfMagErr_g, p.psfMagErr_r, p.psfMagErr_i "
        "FROM PhotoPrimary AS p "
        f"WHERE p.ra BETWEEN {g.ra - dra:.5f} AND {g.ra + dra:.5f} "
        f"AND p.dec BETWEEN {g.dec - half:.5f} AND {g.dec + half:.5f} "
        "AND p.type = 6 AND p.clean = 1 "
        "AND p.psfMag_r BETWEEN 14.5 AND 20.0"
    )
    return SKYSERVER_SQL + "?" + urlencode({"cmd": sql, "format": "csv"})


def _gaia_query(g: Galaxy) -> str:
    radius = g.fov_arcmin / 60.0 * 0.72  # circumscribes the square cutout
    adql = (
        "SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.bp_rp, g.parallax, "
        "g.parallax_over_error, g.pm, g.pmra, g.pmdec, g.pmra_error, g.pmdec_error, g.ruwe, "
        "s.g_sdss_mag, s.r_sdss_mag, s.i_sdss_mag, s.g_sdss_flag, s.r_sdss_flag, s.i_sdss_flag "
        "FROM gaiadr3.gaia_source AS g "
        "LEFT OUTER JOIN gaiadr3.synthetic_photometry_gspc AS s ON s.source_id = g.source_id "
        f"WHERE 1 = CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', {g.ra:.5f}, {g.dec:.5f}, {radius:.4f}))"
    )
    return GAIA_TAP_SYNC + "?" + urlencode({"REQUEST": "doQuery", "LANG": "ADQL",
                                           "FORMAT": "csv", "QUERY": adql})


def _rc3_query(g: Galaxy) -> str:
    """Every RC3 galaxy in the frame: the target plus any companions."""
    return VIZIER_ASU + "?" + urlencode({
        "-source": "VII/155/rc3",
        "-out": "name,altname,PGC,RA2000,DE2000,D25,R25,PA,BT,B-VT,U-BT,T",
        "-out.max": "200",
        "-c": f"{g.ra:.5f}{g.dec:+.5f}",
        "-c.rs": f"{g.fov_arcmin * 60 * 0.72:.0f}",
    })


SVO_FILTER = "http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id="
CALSPEC = "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/"
PICKLES = "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/pickles/dat_uvk/"
CIE_FILES = "https://files.cie.co.at/Publications-datasets/"

#: SVO Filter Profile Service identifiers. SDSS curves are Doi et al. (2010):
#: the response of the whole system -- atmosphere at airmass 1.3, telescope,
#: filter and CCD -- as used by the survey's own calibration.
FILTER_IDS = {
    "sdss_u": "SLOAN/SDSS.u", "sdss_g": "SLOAN/SDSS.g", "sdss_r": "SLOAN/SDSS.r",
    "sdss_i": "SLOAN/SDSS.i", "sdss_z": "SLOAN/SDSS.z",
    "galex_fuv": "GALEX/GALEX.FUV", "galex_nuv": "GALEX/GALEX.NUV",
    "acs_f435w": "HST/ACS_WFC.F435W", "acs_f555w": "HST/ACS_WFC.F555W",
    "acs_f658n": "HST/ACS_WFC.F658N", "acs_f814w": "HST/ACS_WFC.F814W",
    "irac_1": "Spitzer/IRAC.I1", "irac_4": "Spitzer/IRAC.I4", "mips_24": "Spitzer/MIPS.24mu",
}

#: Pickles (1998) stellar spectral library: the colour of each kind of star,
#: to read galaxy colours against. Index -> spectral type.
PICKLES_STARS = {9: "A0V", 12: "A5V", 14: "F0V", 16: "F5V", 23: "G0V", 26: "G2V",
                 31: "K0V", 36: "K5V", 38: "M0V", 44: "M5V", 3: "B0V", 6: "B5-7V",
                 78: "K0III", 95: "M0III"}

SINGS = "https://irsa.ipac.caltech.edu/data/SPITZER/SINGS/galaxies/ngc5194/"
THINGS = "https://things.www3.mpia.de/Data_files/"
HST_M51 = "https://archive.stsci.edu/pub/hlsp/m51/version1/"


def _build_sources() -> list[Source]:
    src: list[Source] = []
    # --- colour science reference data
    for name, fid in FILTER_IDS.items():
        src.append(Source(f"reference/filter_{name}.dat", SVO_FILTER + fid,
                          f"Filter transmission curve {fid}",
                          "SVO Filter Profile Service (Rodrigo et al. 2012, 2020)",
                          kind="text"))
    src += [
        Source("reference/CIE_xyz_1931_2deg.csv", CIE_FILES + "CIE_xyz_1931_2deg.csv",
               "CIE 1931 colour-matching functions (2-degree observer), 360-830 nm",
               "CIE 2019, Colour-matching functions of CIE 1931 standard colorimetric "
               "observer; DOI 10.25039/CIE.DS.xvudnb9b", kind="text"),
        Source("reference/CIE_std_illum_D65.csv", CIE_FILES + "CIE_std_illum_D65.csv",
               "CIE standard illuminant D65 (the white point of sRGB displays)",
               "CIE 2019, CIE standard illuminant D65; DOI 10.25039/CIE.DS.hjfjmt59", kind="text"),
        Source("reference/CIE_sle_scotopic.csv", CIE_FILES + "CIE_sle_scotopic.csv",
               "CIE 1951 scotopic luminous efficiency V'(lambda): what rods respond to",
               "CIE 2019, CIE spectral luminous efficiency for scotopic vision; "
               "DOI 10.25039/CIE.DS.gr6w4b5g", kind="text"),
        Source("reference/sun_reference_stis_002.fits", CALSPEC + "sun_reference_stis_002.fits",
               "The Sun's measured spectrum, 0.12-2.7 microns (HST/CALSPEC composite)",
               "Bohlin, Gordon & Tremblay 2014, PASP 126, 711 (CALSPEC)"),
        Source("reference/alpha_lyr_stis_012.fits", CALSPEC + "alpha_lyr_stis_012.fits",
               "Vega's measured spectrum (HST/CALSPEC)",
               "Bohlin, Gordon & Tremblay 2014, PASP 126, 711 (CALSPEC)"),
        Source("reference/brown2014_galsedatlas.fits",
               "https://archive.stsci.edu/hlsps/galsedatlas/hlsp_galsedatlas_multi_multi_all_multi_v1_cat.fits",
               "Brown et al. atlas: integrated spectra of 129 nearby galaxies (optical part from "
               "drift-scan spectrophotometry), plus their matched-aperture SDSS photometry",
               "Brown et al. 2014, ApJS 212, 18; MAST HLSP GALSEDATLAS, DOI 10.17909/t9-5bxk-dh29"),
    ]
    for idx, sp in PICKLES_STARS.items():
        src.append(Source(f"reference/pickles_uk_{idx}.fits", PICKLES + f"pickles_uk_{idx}.fits",
                          f"Pickles library spectrum of a {sp} star",
                          "Pickles 1998, PASP 110, 863"))

    # --- catalogues of each field
    for g in GALAXIES:
        src.append(Source(
            f"catalogs/{g.key}_sdss_stars.csv", _sdss_star_query(g),
            f"{g.name}: SDSS DR17 point sources with PSF magnitudes in the field",
            "SDSS DR17 (Abdurro'uf et al. 2022)", input=False, kind="csv", pause=1.0))
        src.append(Source(
            f"catalogs/{g.key}_gaia.csv", _gaia_query(g),
            f"{g.name}: Gaia DR3 sources (foreground-star masks) with GSPC synthetic SDSS "
            "magnitudes",
            "Gaia DR3 (Gaia Collaboration, Vallenari et al. 2023); GSPC (Montegriffo et al. 2023)",
            kind="csv", pause=1.0))
        src.append(Source(
            f"catalogs/{g.key}_dust.xml",
            IRSA_DUST + "?" + urlencode({"locstr": f"{g.ra:.5f} {g.dec:+.5f} equ j2000", "regSize": "2.0"}),
            f"{g.name}: Milky Way dust reddening E(B-V) in front of the galaxy",
            "Schlegel, Finkbeiner & Davis 1998; Schlafly & Finkbeiner 2011 (IRSA DUST service)",
            kind="xml", pause=0.5))
        src.append(Source(
            f"catalogs/{g.key}_rc3.tsv", _rc3_query(g),
            f"{g.name}: RC3 isophotal size (D25), axis ratio and position angle of every "
            "catalogued galaxy in the frame",
            "de Vaucouleurs et al. 1991 (RC3), via VizieR VII/155", kind="tsv", pause=0.5))

    # --- the gallery: calibrated images (in tiles), and the survey's own
    # rendering of one tile, to check we read the pixels as the survey does
    for g in GALAXIES:
        for t in tile_layout(g):
            src.append(Source(
                t.filename, tile_url(t, g.pixscale),
                f"{g.name}: SDSS g, r, i images, tile ({t.iy}, {t.ix}) of the "
                f"{g.fov_arcmin:g}' frame at {g.pixscale}\"/pixel (nanomaggies per native pixel)",
                "SDSS DR8 imaging (Aihara et al. 2011) via the Legacy Surveys viewer (Dey et al. 2019)",
                pause=1.0))
        c = central_tile(g)
        src.append(Source(
            f"cutouts/{g.key}/tile_{c.iy}_{c.ix}_viewer.jpg", tile_url(c, g.pixscale, kind="jpeg"),
            f"{g.name}: the viewer's own colour JPEG of the central tile's pixels",
            "Legacy Surveys viewer (Dey et al. 2019)", input=False, kind="jpeg", pause=1.0))
    for key in ("m81", "m82"):
        g = galaxy_by_key(key)
        src.append(Source(
            f"cutouts/{g.key}_sdss_single.fits", cutout_url(g),
            f"{g.name}: the whole frame requested as one cutout, to check the tile stitching",
            "SDSS DR8 imaging via the Legacy Surveys viewer", input=False, pause=1.0))

    # --- M51 across the spectrum
    m51 = galaxy_by_key("m51")
    src += [
        Source("multiwavelength/m51_galex.fits", cutout_url(m51, layer="galex", pixscale=1.4, size=800),
               "M51: GALEX near- and far-ultraviolet images on the optical grid",
               "GALEX (Martin et al. 2005) via the Legacy Surveys viewer", pause=1.0),
        Source("multiwavelength/m51_lsdr9north.fits",
               cutout_url(m51, layer="ls-dr9", pixscale=1.4, size=686),
               "M51 in Legacy Surveys DR9-north g, r, z: a pipeline cautionary tale",
               "Legacy Surveys DR9 (Dey et al. 2019; Zou et al. 2019)", input=False, pause=1.0),
        Source("multiwavelength/m51_lsdr9north_viewer.jpg",
               cutout_url(m51, layer="ls-dr9", pixscale=1.4, size=686, kind="jpeg"),
               "The viewer's own rendering of the DR9-north M51 cutout",
               "Legacy Surveys viewer", input=False, kind="jpeg", pause=1.0),
        Source("multiwavelength/ngc5194_irac1.fits", SINGS + "IRAC/ngc5194_v7.phot.1.fits",
               "M51: Spitzer IRAC 3.6 micron mosaic (old stars)",
               "SINGS DR5 (Kennicutt et al. 2003; Dale et al. 2007)"),
        Source("multiwavelength/ngc5194_irac4.fits", SINGS + "IRAC/ngc5194_v7.phot.4.fits",
               "M51: Spitzer IRAC 8.0 micron mosaic (glowing dust: PAH molecules)",
               "SINGS DR5 (Kennicutt et al. 2003; Dale et al. 2007)"),
        Source("multiwavelength/ngc5194_mips24.fits", SINGS + "MIPS/ngc5194_mips24_crop_v5-0.fits",
               "M51: Spitzer MIPS 24 micron image (warm dust heated by young stars)",
               "SINGS DR5 (Kennicutt et al. 2003; Dale et al. 2007)"),
        Source("multiwavelength/ngc5194_halpha_sub.fits", SINGS + "Ancillary/ngc5194_HA_SUB_dr4.fits",
               "M51: continuum-subtracted H-alpha image (ionised gas around the youngest stars)",
               "SINGS ancillary data, KPNO 2.1m (Kennicutt et al. 2003)"),
        Source("multiwavelength/ngc5194_things_mom0.fits", THINGS + "NGC_5194_NA_MOM0_THINGS.FITS",
               "M51: VLA 21-cm map of neutral hydrogen (THINGS, natural weighting, moment 0)",
               "Walter et al. 2008, AJ 136, 2563 (THINGS)"),
    ]
    for band, filt in (("b", "F435W"), ("v", "F555W"), ("i", "F814W"), ("h", "F658N")):
        src.append(Source(f"multiwavelength/h_m51_{band}_s20_drz_sci.fits",
                          HST_M51 + f"h_m51_{band}_s20_drz_sci.fits",
                          f"M51: Hubble ACS mosaic in {filt}, block-averaged to 0.2\"/pixel",
                          "Mutchler et al. 2005 (Hubble Heritage M51 mosaic, MAST HLSP)"))
    return src


SOURCES: list[Source] = _build_sources()


# --------------------------------------------------------------------------
# Published values we grade ourselves against (never used as inputs)
# --------------------------------------------------------------------------

PUBLISHED = {
    "sun_abs_mag_ab": {
        "reference": "Willmer 2018, ApJS 236, 47 (Table 3)",
        "values": {"u": 6.39, "g": 5.11, "r": 4.65, "i": 4.53, "z": 4.50},
    },
}

#: The rectangular apertures Brown et al. (2014) measured their SDSS
#: photometry in (their Table 2): (length along the position angle, width,
#: position angle in degrees east of north), in arcsec. They trace the regions
#: covered by the drift-scan spectra -- mostly 55"-wide strips along the major
#: axis -- so they are not whole-galaxy apertures. Centres are the positions
#: in their Table 1. Re-measuring our images in exactly these rectangles
#: compares like with like.
BROWN_APERTURES = {
    "NGC_0628": (346.0, 55.0, 70.0),
    "NGC_4486": (300.0, 300.0, 90.0),
    "NGC_4594": (287.0, 55.0, 113.0),
    "NGC_4631": (512.0, 56.0, 101.0),
    "NGC_4676_A": (40.0, 155.0, 90.0),
    "NGC_4826": (330.0, 55.0, 120.0),
    "NGC_5055": (416.0, 56.0, 120.0),
    "NGC_5194": (370.0, 56.0, 158.0),
    "NGC_5195": (190.0, 56.0, 118.0),
    "NGC_5866": (154.0, 55.0, 120.0),
}
