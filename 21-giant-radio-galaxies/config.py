"""Project configuration: paths, data sources, the giants, and every assumption.

Every number this project *chooses* rather than measures lives here: the
cosmology that turns angles into megaparsecs, the contour level that decides
where a radio galaxy ends, the threshold that makes a radio galaxy a "giant",
the priors of the selection model, the physics constants of synchrotron
ageing. Each sits next to the reason for it and, where it comes from a paper,
the paper.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlencode

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
CATALOG_DIR = DATA_DIR / "catalogs"
IMAGE_DIR = DATA_DIR / "images"
SDSS_DIR = DATA_DIR / "sdss"
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, CATALOG_DIR, IMAGE_DIR, SDSS_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-project-atlas/I18 (hobby research; python-urllib)"


# --------------------------------------------------------------------------
# Cosmology: how an angle on the sky becomes a length
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Cosmology:
    """A flat Lambda-CDM universe: everything an angle-to-length conversion needs."""

    name: str
    H0: float  # km/s/Mpc
    Om0: float  # matter density today; dark energy is 1 - Om0 (flat, radiation ignored)
    reference: str


#: The parameters every giant-radio-galaxy catalogue used here adopts (Oei et
#: al. 2023, 2024; Mostert et al. 2024): Planck Collaboration 2020 (A&A 641,
#: A6), Table 2, TT,TE,EE+lowE+lensing+BAO: h = 0.6766, Omega_M = 0.3111, flat,
#: no radiation term. With them we reproduce both catalogues' megaparsec
#: lengths from their angles and redshifts (checked in ``lengths.py``).
PLANCK18 = Cosmology("Planck 2018", 67.66, 0.3111,
                     "Planck Collaboration 2020, A&A 641, A6 (Table 2, TT,TE,EE+lowE+lensing+BAO)")

#: This atlas's own cosmology from the cosmic microwave background (project D3,
#: ``09-planck-cmb``): a CAMB fit to the Planck 2018 TT spectrum gave
#: H0 = 66.66 +- 0.90 km/s/Mpc and Omega_m = 0.3247 +- 0.0132.
D3_CMB = Cosmology("D3 (our CMB fit)", 66.66, 0.3247,
                   "this atlas, D3: CAMB fit to the Planck 2018 TT spectrum")

#: This atlas's supernova cosmology (project D1, ``05-supernova-cosmology``):
#: Cepheid-calibrated H0 = 73.55 +- 1.02 and flat Omega_m = 0.332 +- 0.018.
#: The two sides of the Hubble tension change every length by ~10%.
D1_SNE = Cosmology("D1 (our supernova fit)", 73.55, 0.332,
                   "this atlas, D1: Pantheon+ with SH0ES Cepheid calibration")

COSMOLOGY = PLANCK18
COSMOLOGIES = (PLANCK18, D3_CMB, D1_SNE)

#: A radio galaxy is a "giant" when its projected proper length exceeds this.
#: 0.7 Mpc is the modern convention (it was 1 Mpc for H0 = 50 in the 1970s,
#: which is the same physical size); all catalogues used here adopt it.
GIANT_THRESHOLD_MPC = 0.7


# --------------------------------------------------------------------------
# The length distribution: Oei et al. (2023) forward model
# --------------------------------------------------------------------------

#: Oei et al. (2023, A&A 672, A163) found 2,060 giants in LoTSS DR2 by eye and
#: inferred the tail index of the intrinsic length distribution from those at
#: z < 0.5 with a definite host (1,473 giants; 811 at z < 0.25). We redo that
#: inference on the same giants with our own implementation of their model:
#: Pareto-distributed lengths, random orientations, a "fuzzy" angular-size
#: threshold around the 5' their search used, and lobe surface-brightness
#: selection. Settings below are theirs (their Sects. 2.8 and 4.1).
LENGTH_MODEL = {
    "z_max": (0.5, 0.25),
    "phi_threshold_arcmin": 5.0,  # their search's intended angular-length threshold
    "l_ref_mpc": 0.7,  # reference length for the surface-brightness model
    "nu_obs_mhz": 144.0,
    "b_threshold_jy_deg2": 25.0,  # 1 sigma of LoTSS DR2 at 6" (Shimwell et al. 2022)
    "lobe_spectral_index": -1.0,
    "zeta": -2.0,  # self-similar growth: surface brightness ~ length^-2
    # Gaussian priors from their own measurements (their Sect. 4.1.3):
    "b_ref_prior_jy_deg2": (1000.0, 750.0),
    "sigma_ref_prior": (1.3, 0.65),
    # Published posterior (their Table 4) to compare with.
    "published": {
        0.5: {"xi": (-3.5, 0.5), "xi_map": -3.55, "half_width_arcmin": (1.9, 0.2),
              "b_ref": (720.0, 380.0), "sigma_ref": (1.2, 0.2), "n": 1473},
        0.25: {"xi": (-3.5, 0.4), "xi_map": -3.5, "half_width_arcmin": (1.7, 0.3),
               "b_ref": (1020.0, 490.0), "sigma_ref": (1.3, 0.4), "n": 811},
    },
}

#: Mostert et al. (2024, A&A 691, A185, their Eqs. 7-9) let the power law bend:
#: f(l_p) ~ (l_p / 0.7 Mpc)^xi(l_p) with xi(l_p) = xi_0 + dxi (l_p - 0.7)/(5 - 0.7),
#: and found xi_0 = -2.8 +- 0.2, dxi = -2.4 +- 0.3 from 5,000+ giants at z < 0.5.
MOSTERT_LENGTH_MODEL = {"xi_at_threshold": (-2.8, 0.2), "delta_xi": (-2.4, 0.3), "l_p2_mpc": 5.0}


# --------------------------------------------------------------------------
# Radio images
# --------------------------------------------------------------------------

#: LoTSS cutout API: cuts from the survey's own mosaics, in Jy/beam. 6"
#: resolution with 1.5" pixels, or the 20" "low-resolution" images (4.5"
#: pixels) that show faint, extended lobes better. The service refuses cutouts
#: larger than 2 degrees on a side. DR2 (Shimwell et al. 2022, A&A 659, A1)
#: covers 27% of the northern sky and is what every giant catalogue was built
#: from; DR3 (Shimwell et al. 2026, arXiv:2602.15949) covers 88% of it, with
#: 9" instead of 6" resolution south of declination +10 deg.
LOTSS_URL = {"dr2": {6: "https://lofar-surveys.org/dr2-cutout.fits",
                     20: "https://lofar-surveys.org/dr2-low-cutout.fits"},
             "dr3": {6: "https://lofar-surveys.org/dr3-cutout.fits",
                     20: "https://lofar-surveys.org/dr3-low-cutout.fits"}}
LOTSS_MAX_SIZE_ARCMIN = 120.0
#: Giants larger than one cutout are fetched as a 2 x 2 mosaic of tiles this
#: size, offset by +-0.36 tile from the host, and stitched on our own grid.
MOSAIC_TILE_ARCMIN = 100.0
LOTSS_FREQ_MHZ = 144.0
#: LoTSS DR2's flux scale is tied to Roger et al. (1973) / Scaife & Heald (2012)
#: to about 10% (Shimwell et al. 2022). NVSS is good to ~3% (Condon et al. 1998).
LOTSS_FLUX_SCALE_ERR = 0.10
NVSS_FLUX_SCALE_ERR = 0.03

#: NVSS at 1.4 GHz, 45" resolution (Condon et al. 1998, AJ 115, 1693). NRAO's
#: own postage-stamp server answers with HTTP 500 (October 2026), so we take
#: NVSS from NASA's SkyView, which resamples the survey's own images in Jy/beam.
#: NVSS observed in snapshots with the VLA's most compact configuration: it
#: cannot see smooth structure larger than about 15' (Condon et al. 1998), so
#: its fluxes of the biggest lobes are lower limits.
SKYVIEW_URL = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
NVSS_FREQ_MHZ = 1400.0
NVSS_BEAM_ARCSEC = 45.0
NVSS_LARGEST_SCALE_ARCMIN = 15.0

#: Legacy Surveys sky viewer (Dey et al. 2019, AJ 157, 168): colour JPEGs of
#: DR10 g, r, z imaging (DR9 in the north, where DR10 adds nothing but colour
#: artefacts) for the optical backdrop and the host-galaxy close-ups.
LEGACY_URL = "https://www.legacysurvey.org/viewer/jpeg-cutout"
LEGACY_LAYER = "ls-dr10"
LEGACY_MAX_PX = 1600
#: Where the Legacy Surveys have no imaging (near the celestial pole), the
#: Pan-STARRS1 colour survey (Chambers et al. 2016) through CDS's hips2fits.
HIPS2FITS_URL = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
PS1_HIPS = "CDS/P/PanSTARRS/DR1/color-z-zg-g"

#: Where a radio galaxy "ends": the 3-sigma contour of the image, as in the
#: catalogues' own images (contours at 3, 5, 10 sigma). Islands smaller than one
#: beam are noise spikes, not lobes.
LAS_SIGMA = 3.0
CONTOUR_SIGMAS = (3, 5, 10, 20, 40, 80, 160)
MIN_ISLAND_BEAMS = 1.0

#: Oei et al. (2023) estimate their angular lengths good to 0.15'.
PUBLISHED_ANGLE_ERR_ARCMIN = 0.15


# --------------------------------------------------------------------------
# The giants we image
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Giant:
    """One giant radio galaxy we fetch images of and measure ourselves."""

    key: str
    name: str
    #: Host galaxy position (deg, J2000) and its name.
    ra: float
    dec: float
    host: str
    #: Redshift we adopt, its uncertainty, kind ("s" spectroscopic, "p"
    #: photometric) and where it comes from.
    z: float
    z_err: float
    z_kind: str
    z_ref: str
    #: Why this one is in the project.
    why: str
    #: Field of view to fetch (arcmin, square). Wide enough for the whole giant
    #: plus empty sky around it to measure the noise.
    fov_arcmin: float
    #: Lobe "seeds": points (deg) inside each end of the radio galaxy, placed by
    #: eye on the LoTSS image -- the step that needs a human, as in every giant
    #: catalogue. The length itself is then measured from the 3-sigma contour.
    seeds: tuple = ()
    #: Extra points (deg) inside emission that belongs to the giant but is not
    #: connected to either end (a core, a jet knot, a detached lobe piece).
    extra_seeds: tuple = ()
    #: Sky positions (deg) of unrelated sources to exclude even if their
    #: emission touches the giant's.
    exclude: tuple = ()
    #: For disputed giants: seeds for a conservative association -- only the
    #: emission clearly tied to the host (its own double or arm) -- to set
    #: against the claimed one.
    alt_seeds: tuple = ()
    #: Resolution (arcsec) at which we measure the length: 20" for giants with
    #: compact, bright lobes; 60" (the 20" image smoothed) for faint, diffuse
    #: lobes that break up into noise at 20".
    measure_beam: float = 20.0
    #: Which group(s) the giant appears in: "gallery" (the keeper picture),
    #: "records" (the claimed longest giants), "ageing" (spectral ages),
    #: "dr3" (giants outside the public DR2 footprint, imaged in DR3).
    groups: tuple = ("gallery",)
    #: Lengths claimed in the literature: (label, angle arcmin or None, Mpc).
    claims: tuple = ()
    #: LoTSS release the giant is imaged in, and where its optical image comes
    #: from ("legacy" or, where the Legacy Surveys have no data, "ps1").
    survey: str = "dr2"
    optical: str = "legacy"
    #: Tiles per side: 1 for a single cutout, 2 for a 2 x 2 mosaic.
    tiles: int = 1

    @property
    def size_arcmin(self) -> float:
        if self.tiles > 1:
            return self.fov_arcmin
        return min(self.fov_arcmin, LOTSS_MAX_SIZE_ARCMIN)


GIANTS: list[Giant] = [
    Giant("porphyrion", "Porphyrion", 232.38410, 60.25960, "SDSS J152932.16+601534.4",
          0.896, 0.001, "s", "Keck/LRIS spectrum (Oei et al. 2024, Nature 633, 537)",
          "the longest jets known: 6.43 Mpc projected, ~7 Mpc in 3D, seen 7.5 Gyr ago",
          22.0, seeds=((232.466, 60.360), (232.300, 60.195)), extra_seeds=((232.3841, 60.2596),),
          groups=("gallery", "records"),
          claims=(("Oei et al. 2024 (Nature)", None, 6.43),
                  ("Andernach et al. 2025", 13.2, 6.16),
                  ("Mostert et al. 2024 (RGZ, z = 0.916)", 17.36, 8.39),
                  ("Hardcastle et al. 2023 (host and z uncertain)", None, 6.9))),
    Giant("alcyoneus", "Alcyoneus", 123.590372, 52.402795, "SDSS J081421.68+522410.0",
          0.24674, 0.00006, "s", "SDSS DR12 (via Oei et al. 2022)",
          "the record holder of 2022: 4.99 Mpc, from an ordinary elliptical galaxy",
          30.0, seeds=((123.430, 52.475), (123.720, 52.315)), extra_seeds=((123.5904, 52.4028),),
          measure_beam=60.0, groups=("gallery", "records", "ageing"),
          claims=(("Oei et al. 2022/2023", 20.8, 4.99), ("Andernach et al. 2025", 20.4, 4.74))),
    Giant("3c236", "3C 236", 151.507278, 34.902906, "2MASX J10060176+3454104",
          0.09936, 0.00002, "s", "SDSS (via Mostert et al. 2024)",
          "the first giant ever found (Willis et al. 1974), 4.3-4.9 Mpc",
          62.0, seeds=((151.300, 35.010), (151.920, 34.695)), extra_seeds=((151.5073, 34.9029),),
          groups=("gallery", "records", "ageing"),
          claims=(("Mostert et al. 2024 (RGZ)", 42.97, 4.88), ("Andernach et al. 2025", 39.0, 4.33),
                  ("Oei et al. 2024 (Nature)", None, 4.6))),
    Giant("j1547", "LEDA 56028", 236.788416, 35.646139, "SDSS J154709.22+353846.1",
          0.0794, 0.00001, "s", "SDSS DR12 (via Oei et al. 2023)",
          "4.08 Mpc and 44' across: a giant larger than the full Moon",
          62.0, seeds=((236.530, 35.775), (237.310, 35.535)), extra_seeds=((236.7884, 35.6461),),
          measure_beam=60.0, groups=("gallery", "ageing"),
          claims=(("Oei et al. 2023", 43.8, 4.08),)),
    Giant("ngc6185", "NGC 6185", 248.3243, 35.3423, "NGC 6185",
          0.03430, 0.00007, "s", "Oei et al. 2023",
          "the largest radio galaxy known to come from a spiral galaxy",
          86.0, seeds=((248.720, 35.300), (247.970, 35.330)), extra_seeds=((248.3243, 35.3423),),
          measure_beam=60.0, claims=(("Oei et al. 2023", 60.0, 2.54),)),
    Giant("ngc315", "NGC 315", 14.4537, 30.3524, "NGC 315",
          0.0165, 0.0001, "s", "Dabhade et al. 2020 (via Mostert et al. 2024)",
          "a nearby giant whose twisting jets span almost a degree of sky",
          110.0, seeds=((15.000, 30.080), (14.170, 30.540)),
          extra_seeds=((14.350, 30.150), (14.4537, 30.3524)), groups=("gallery", "ageing"),
          claims=(("Dabhade et al. 2020", 57.0, 1.19),)),
    Giant("3c31", "3C 31", 16.8540, 32.4125, "NGC 383",
          0.0167, 0.0001, "s", "Dabhade et al. 2020 (via Mostert et al. 2024)",
          "a classic radio galaxy that LOFAR revealed to be a giant (Heesen et al. 2018: 1.1 Mpc)",
          66.0, seeds=((16.630, 32.750), (16.690, 31.990)), extra_seeds=((16.8540, 32.4125),),
          claims=(("Dabhade et al. 2020", 45.0, 0.95), ("Heesen et al. 2018 (LOFAR)", None, 1.1))),
    Giant("j1318", "J1318+2626", 199.597594, 26.439694, "SDSS J131823.42+262622.8",
          0.623, 0.00005, "s", "SDSS DR12 (via Oei et al. 2023)",
          "a giant launched from a quasar: the black hole is also shining brightly",
          18.0, seeds=((199.672, 26.452), (199.495, 26.418)), extra_seeds=((199.5976, 26.4397),),
          groups=("gallery", "records"),
          claims=(("Oei et al. 2023", 11.0, 4.62), ("Andernach et al. 2025", 10.72, 4.37))),
    Giant("j0819", "J0819+3235", 124.98508, 32.593796, "SDSS J081956.41+323537.6",
          0.7535, 0.0005, "s", "Andernach et al. 2025 (spectroscopic)",
          "the longest of Oei et al.'s 2,060 discoveries", 18.0,
          seeds=((125.052, 32.644), (124.952, 32.545)), extra_seeds=((124.98508, 32.593796),),
          groups=("gallery", "records"),
          claims=(("Oei et al. 2023 (photo-z 0.749)", 11.2, 5.07), ("Andernach et al. 2025", 10.56, 4.66),
                  ("Mostert et al. 2024 (RGZ, z = 0.793)", 10.66, 4.93))),
    Giant("ngc2789", "NGC 2789", 138.7486, 29.7303, "NGC 2789",
          0.0212, 0.0001, "s", "Oei et al. 2023",
          "a giant from a lenticular galaxy, 35' across", 50.0,
          seeds=((138.980, 29.720), (138.650, 29.750)), extra_seeds=((138.7486, 29.7303),), measure_beam=60.0,
          claims=(("Oei et al. 2023", 34.5, 0.92),)),
    Giant("j0838", "J0838+5327", 129.5074, 53.4540, "SDSS J083801.77+532714.4",
          0.5445, 0.0001, "s", "Andernach et al. 2025",
          "a claimed rival to Porphyrion (6.58 Mpc), resting on a disputed northern lobe",
          26.0, seeds=((129.598, 53.600), (129.445, 53.350)), extra_seeds=((129.5074, 53.4540),),
          alt_seeds=((129.5074, 53.4540), (129.445, 53.350)), groups=("records",), claims=(("Andernach et al. 2025", 17.2, 6.58),)),
    Giant("iltj0036", "ILTJ003621.03+253500.9", 9.099713, 25.603076, "(RGZ host)",
          0.62542, 0.0, "p", "Mostert et al. 2024", "the second-longest entry in Mostert et al.'s catalogue",
          28.0, seeds=((9.280, 25.690), (8.990, 25.595)), extra_seeds=((9.0997, 25.6031), (9.070, 25.565)),
          alt_seeds=((9.0997, 25.6031), (9.070, 25.565)), measure_beam=60.0, groups=("records",), claims=(("Mostert et al. 2024 (RGZ)", 18.97, 7.98),)),
    Giant("iltj1327", "ILTJ132739.02+350643.0", 201.8999692, 35.1119627, "(RGZ host)",
          0.891, 0.0, "p", "Mostert et al. 2024", "third-longest entry in Mostert et al.'s catalogue",
          20.0, seeds=((201.945, 35.112), (201.735, 35.105)), extra_seeds=((201.860, 35.112),),
          alt_seeds=((201.945, 35.112), (201.860, 35.112)), measure_beam=60.0, groups=("records",), claims=(("Mostert et al. 2024 (RGZ)", 13.04, 6.25),)),
    Giant("iltj0927", "ILTJ092754.27+351037.1", 141.957353, 35.181085, "(RGZ host)",
          0.748, 0.0, "p", "Mostert et al. 2024", "fourth-longest entry in Mostert et al.'s catalogue",
          20.0, seeds=((142.020, 35.168), (141.780, 35.252)), extra_seeds=((141.900, 35.190),),
          alt_seeds=((142.020, 35.168), (141.900, 35.190)), groups=("records",), claims=(("Mostert et al. 2024 (RGZ)", 12.79, 5.79),)),
    Giant("iltj0834", "ILTJ083456.98+541742.4", 128.72598, 54.292234, "(RGZ host)",
          0.83775, 0.0, "p", "Mostert et al. 2024", "fifth-longest entry in Mostert et al.'s catalogue",
          19.0, seeds=((128.890, 54.297), (128.600, 54.297)),
          extra_seeds=((128.775, 54.297), (128.726, 54.292), (128.665, 54.290)),
          alt_seeds=((128.775, 54.297), (128.600, 54.297)), groups=("records",), claims=(("Mostert et al. 2024 (RGZ)", 12.02, 5.66),)),
    Giant("j0103", "J0103+2305", 15.8522, 23.0841, "SDSS J010324.53+230502.8",
          0.265, 0.0001, "s", "Andernach et al. 2025 (lower limit: host uncertain)",
          "a 24' giant whose host is uncertain", 34.0,
          seeds=((16.000, 23.065), (15.615, 23.135)), extra_seeds=((15.8522, 23.0841),),
          groups=("records",),
          claims=(("Andernach et al. 2025", 24.1, 5.9),)),
    Giant("j1108", "J1108+2916", 167.1436, 29.2753, "SDSS J110834.45+291630.9",
          0.80, 0.1, "p", "Andernach et al. 2025 (photometric)", "a straight 11.6' giant",
          18.0, seeds=((167.190, 29.337), (167.075, 29.180)),
          extra_seeds=((167.130, 29.290), (167.100, 29.245), (167.1436, 29.2753)),
          groups=("records",), claims=(("Andernach et al. 2025", 11.56, 5.2),)),
]


#: Famous giants outside the public DR2 footprint, imaged in LoTSS DR3. Hosts,
#: redshifts and published sizes from the catalogue of Mostert et al. (2024),
#: which credits the original discoveries.
DR3_GIANTS: list[Giant] = [
    Giant("ngc2300", "NGC 2300", 113.085358, 85.708862, "NGC 2300",
          0.00635, 0.0001, "s", "Oei et al. 2023",
          "the largest radio galaxy on the sky: 2.2 degrees, four full Moons end to end", 170.0,
          seeds=((114.2398, 85.3417), (101.8375, 86.4737)), extra_seeds=((112.6279, 85.7687), (107.0058, 86.2088)),
          groups=("dr3",), survey="dr3", optical="ps1", tiles=2, measure_beam=90.0,
          claims=(("Oei et al. 2023", 132.0, 1.07),)),
    Giant("ngc6251", "NGC 6251", 248.1332, 82.5379, "NGC 6251",
          0.0240, 0.0001, "s", "Dabhade et al. 2020 (via Mostert et al. 2024)",
          "a classic giant with one of the straightest, brightest jets in the sky", 95.0,
          seeds=((245.1823, 82.5297), (251.0405, 82.2848)), extra_seeds=((248.1332, 82.5379),),
          groups=("dr3",), survey="dr3", claims=(("Dabhade et al. 2020", 65.7, 1.97),)),
    Giant("da240", "DA 240", 117.15375, 55.816111, "DA 240 host (4C 56.16)",
          0.0356, 0.0001, "s", "Dabhade et al. 2020 (via Mostert et al. 2024)",
          "one of the first giants known (Willis et al. 1974)", 52.0,
          seeds=((117.5792, 55.9729), (116.7046, 55.6927)), extra_seeds=((117.15375, 55.816111),),
          groups=("dr3",), survey="dr3", claims=(("Dabhade et al. 2020", 35.15, 1.54),)),
    Giant("4c7308", "4C 73.08", 147.4415, 73.2398, "4C 73.08 host",
          0.0581, 0.0001, "s", "Dabhade et al. 2020 (via Mostert et al. 2024)",
          "a classic giant with lobes at both ends", 28.0,
          seeds=((147.9603, 73.2914), (147.0441, 73.2594)), extra_seeds=((147.4415, 73.2398),),
          groups=("dr3",), survey="dr3", claims=(("Dabhade et al. 2020", 16.0, 1.11),)),
    Giant("3c326", "3C 326", 238.038333, 20.089722, "3C 326 N",
          0.0895, 0.0001, "s", "Dabhade et al. 2020 (via Mostert et al. 2024)",
          "a 2-Mpc giant from a galaxy with almost no star formation", 32.0,
          seeds=((238.1572, 20.1177), (237.8288, 20.0571)), extra_seeds=((238.038333, 20.089722),),
          groups=("dr3",), survey="dr3", claims=(("Dabhade et al. 2020", 20.1, 2.08),)),
    Giant("j1526", "J1526+2620", 231.644896, 26.334223, "2MASX J15263474+2620026",
          0.15072, 0.00002, "s", "SDSS DR12 (via Oei et al. 2023)",
          "the longest giant of Oei et al.'s search that DR2's public images miss: 4.56 Mpc claimed", 42.0,
          seeds=((231.7634, 26.4937), (231.6221, 26.3148)), measure_beam=120.0,
          groups=("dr3",), survey="dr3",
          claims=(("Oei et al. 2023", 28.0, 4.56), ("Andernach et al. 2025", 28.0, 4.40))),
    Giant("j1129", "J1129+2733", 172.300621, 27.553881, "SDSS J112912.14+273313.9",
          0.07324, 0.00001, "s", "SDSS DR12 (via Oei et al. 2023)",
          "a 3.3-Mpc giant almost 40' across", 58.0,
          seeds=((172.5419, 27.7190), (172.0696, 27.2620)), extra_seeds=((172.300621, 27.553881),),
          measure_beam=60.0, groups=("dr3",), survey="dr3",
          claims=(("Oei et al. 2023", 38.6, 3.34), ("Andernach et al. 2025", 40.0, 3.34))),
    Giant("ngc1044", "NGC 1044", 40.2758, 8.7380, "NGC 1044",
          0.0212, 0.0001, "s", "Oei et al. 2023",
          "a Moon-sized giant south of the DR2 footprint (DR3 resolution 9\" here)", 48.0,
          seeds=((40.3747, 8.9237), (40.2419, 8.4138)), extra_seeds=((40.2758, 8.7380),), measure_beam=60.0,
          groups=("dr3",), survey="dr3", claims=(("Oei et al. 2023", 32.7, 0.87),)),
    Giant("ngc7385", "NGC 7385", 342.4774, 11.6086, "NGC 7385",
          0.0262, 0.0001, "s", "Oei et al. 2023",
          "a Moon-sized giant in a galaxy cluster", 72.0,
          seeds=((342.5738, 11.5429), (342.0779, 11.4402)), extra_seeds=((342.4774, 11.6086),),
          measure_beam=60.0, groups=("dr3",), survey="dr3",
          claims=(("Oei et al. 2023", 32.0, 1.04),)),
    Giant("3c35", "3C 35", 18.0094, 49.4764, "3C 35 host",
          0.0677, 0.0001, "s", "Dabhade et al. 2020 (via Mostert et al. 2024)",
          "a 1-Mpc giant whose jets have switched off", 22.0,
          seeds=((18.0304, 49.5658), (18.0162, 49.3832)), extra_seeds=((18.0094, 49.4764),),
          groups=("dr3",), survey="dr3",
          claims=(("Dabhade et al. 2020", 12.5, 1.00),)),
]
GIANTS = GIANTS + DR3_GIANTS


def giant_by_key(key: str) -> Giant:
    for g in GIANTS:
        if g.key == key:
            return g
    raise KeyError(key)


def giants_in(group: str) -> list[Giant]:
    return [g for g in GIANTS if group in g.groups]


# --------------------------------------------------------------------------
# Sizes for comparison
# --------------------------------------------------------------------------

#: The Milky Way's stellar disc, ~30 kpc across (the "100,000 light-years" of
#: the textbooks; the disc has no sharp edge). Andromeda's disc to its
#: 25 mag/arcsec^2 isophote: RC3 D25 = 190.5' at 761 +- 11 kpc (Li et al. 2021,
#: ApJ 920, 84) = 42 kpc. Their separation, 761 kpc, is the size of the whole
#: Local Group's core.
MILKY_WAY_DIAMETER_KPC = 30.0
ANDROMEDA_DIAMETER_KPC = 42.2
ANDROMEDA_DISTANCE_KPC = 761.0
#: The full Moon spans about 31' (it varies from 29.4' to 34.1').
MOON_DIAMETER_ARCMIN = 31.0


# --------------------------------------------------------------------------
# Black holes, galaxies and their environment
# --------------------------------------------------------------------------

#: Kormendy & Ho (2013, ARA&A 51, 511) Eq. 7, for classical bulges and
#: ellipticals: M_BH / 1e9 Msun = 0.309 (sigma / 200 km/s)^4.38, intrinsic
#: scatter 0.29 dex. Oei et al. (2023) used the same relation.
MSIGMA = {"norm_1e9": 0.309, "slope": 4.38, "scatter_dex": 0.29}
#: SDSS velocity dispersions below ~70 km/s are under the spectrographs'
#: resolution and unreliable; above 420 km/s the SDSS pipeline flags them.
#: We also need a spectrum good enough that the dispersion's error is < 20%.
SIGMA_RANGE_KMS = (70.0, 420.0)
SIGMA_MAX_REL_ERR = 0.20

#: SkyServer (SDSS DR17) for spectra near each host, in batches small enough
#: that the query fits in a URL (~100 positions).
SKYSERVER_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
SKYSERVER_BATCH = 100
SPEC_MATCH_ARCMIN = 0.05  # 3": a host and its fibre

#: Radio galaxies that are not giants, for comparison: LoTSS DR2 value-added
#: catalogue sources (Hardcastle et al. 2023, A&A 678, A151) that are resolved,
#: have a redshift, and are radio-loud enough at 144 MHz that star formation
#: cannot explain them (> 1e24 W/Hz; Oei et al. 2024 used the same cut).
VAC_TABLE = "J/A+A/678/A151/catalog"
RADIO_LOUD_W_HZ = 1e24
VAC_COLUMNS = ("Source", "RAJ2000", "DEJ2000", "Ftot", "e_Ftot", "SCode", "LAS", "r_LAS",
               "CompSize", "CompPA", "Assoc", "RAIdeg", "DEIdeg", "zbest", "r_zbest",
               "zspecSDSS", "zDESI", "Size", "L144", "logMass", "q_logMass",
               "rmagrest", "Kmagrest", "W1mag", "PlateSDSS", "MJDSDSS", "FiberSDSS")

#: The SDSS Main Galaxy Sample (as project D2 used it): galaxies brighter than
#: Petrosian r = 17.77, a flux-limited map of the cosmic web to z ~ 0.2. We
#: count neighbours around each radio galaxy in a cylinder.
MGS_TARGET_BITS = 0x40 | 0x80 | 0x100
MGS_R_LIMIT = 17.77
MGS_Z_RANGE = (0.002, 0.25)
#: Neighbour counting: MGS galaxies in a cylinder of 2 Mpc projected radius
#: and +-1000 km/s. The MGS is flux limited, so counts are only ever compared
#: between objects at the same redshift. Redshift range as Oei et al. (2024):
#: below 0.16, where the MGS still maps the cosmic web densely.
ENV_RADIUS_MPC = 2.0
ENV_DV_KMS = 1000.0
ENV_Z_RANGE = (0.02, 0.16)


# --------------------------------------------------------------------------
# Synchrotron physics for lobe ages
# --------------------------------------------------------------------------

#: Spectral ageing in the Jaffe & Perola (1973) model: electrons injected with a
#: power law (spectral index alpha_inj) lose energy to synchrotron and to
#: inverse-Compton scattering of the CMB; the spectrum steepens above a break
#: frequency that falls with time. Assumptions, all standard and all uncertain:
#: injection index -0.6 (hotspot spectra; e.g. Harwood et al. 2013 find
#: -0.5 to -0.8), equal energy in protons and electrons (k = 1), lobes filled
#: uniformly, field at minimum energy for electrons radiating 10 MHz-100 GHz.
ALPHA_INJ = -0.6
EQUIPARTITION = {"k": 1.0, "filling": 1.0, "nu_min_mhz": 10.0, "nu_max_mhz": 1.0e5}
#: CMB energy density as an equivalent magnetic field: 3.25 (1+z)^2 microgauss.
B_CMB0_UG = 3.25


# --------------------------------------------------------------------------
# Catalogue files (CDS) and other fixed downloads
# --------------------------------------------------------------------------

CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp"
TAP_VIZIER = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"


@dataclass
class Source:
    """One file we download, and why."""

    filename: str  # path under data/
    url: str
    what: str
    reference: str
    input: bool = True  # False: a published result used only to grade our own
    kind: str = "text"  # payload check: "fits", "jpeg", "csv", "text", "gz"
    pause: float = 0.0
    timeout: float = 300.0
    params: dict = field(default_factory=dict)


def _cds(catalog: str, filename: str, local: str, what: str, reference: str, **kw) -> Source:
    return Source(f"catalogs/{local}", f"{CDS_FTP}/{catalog}/{filename}", what, reference, **kw)


CATALOG_SOURCES: list[Source] = [
    _cds("J/A+A/672/A163", "ReadMe", "oei2023_ReadMe", "description of the 2,060-giant table",
         "Oei et al. 2023, A&A 672, A163"),
    _cds("J/A+A/672/A163", "tablef1.dat", "oei2023_tablef1.dat",
         "the 2,060 giants Oei et al. found by eye in LoTSS DR2: host, redshift, angular and projected length",
         "Oei et al. 2023, A&A 672, A163"),
    _cds("J/A+A/691/A185", "ReadMe", "mostert2024_ReadMe", "description of the 11,585-giant catalogue",
         "Mostert et al. 2024, A&A 691, A185"),
    _cds("J/A+A/691/A185", "catalog.dat.gz", "mostert2024_catalog.dat.gz",
         "all 11,585 known giants (LoTSS DR2 crowd-sourced + machine learning + literature)",
         "Mostert et al. 2024, A&A 691, A185", kind="gz"),
    _cds("J/A+A/686/A137", "ReadMe", "oei2024_ReadMe", "description of the Cosmic Web tables",
         "Oei et al. 2024, A&A 686, A137"),
    _cds("J/A+A/686/A137", "tab1grs.dat", "oei2024_tab1grs.dat",
         "BORG SDSS Cosmic Web density and dynamical class around luminous giants (z < 0.16)",
         "Oei et al. 2024, A&A 686, A137"),
    _cds("J/A+A/686/A137", "tab1rg.dat", "oei2024_tab1rg.dat",
         "the same for a control sample of general LoTSS DR1 radio galaxies",
         "Oei et al. 2024, A&A 686, A137"),
    _cds("J/A+A/699/A257", "ReadMe", "andernach2025_ReadMe", "description of the >3 Mpc table",
         "Andernach et al. 2025, A&A 699, A257"),
    _cds("J/A+A/699/A257", "tablea1.dat", "andernach2025_tablea1.dat",
         "142 giants larger than 3 Mpc, each re-examined by eye (hosts, redshifts, angular sizes revised)",
         "Andernach et al. 2025, A&A 699, A257", input=False),
]


def vac_query(dec_lo: float, dec_hi: float) -> str:
    """ADQL for one declination strip of the LoTSS DR2 value-added catalogue:
    resolved sources with a redshift (photometric or spectroscopic)."""
    cols = ", ".join(f'"{c}"' for c in VAC_COLUMNS)
    return (f'SELECT {cols} FROM "{VAC_TABLE}" WHERE "Resolved" = 1 AND "zbest" IS NOT NULL '
            f'AND "DEJ2000" >= {dec_lo} AND "DEJ2000" < {dec_hi}')


#: Declination strips keep each TAP answer small enough to come back in one go.
VAC_DEC_STRIPS = ((-5, 25), (25, 32), (32, 38), (38, 45), (45, 52), (52, 60), (60, 91))


def vac_sources() -> list[Source]:
    out = []
    for lo, hi in VAC_DEC_STRIPS:
        params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": "500000",
                  "QUERY": vac_query(lo, hi)}
        out.append(Source(f"catalogs/vac_resolved_dec{lo:+03d}_{hi:+03d}.csv",
                          f"{TAP_VIZIER}?{urlencode(params)}",
                          "LoTSS DR2 value-added catalogue: resolved sources with redshifts, sizes, luminosities, host masses",
                          "Hardcastle et al. 2023, A&A 678, A151", kind="csv", pause=1.0, timeout=900))
    return out


def lotss_url(g: Giant, resolution: int, *, release: str | None = None, pos: tuple | None = None,
              size_arcmin: float | None = None) -> str:
    ra, dec = pos if pos is not None else (g.ra, g.dec)
    size = size_arcmin if size_arcmin is not None else g.size_arcmin
    return f"{LOTSS_URL[release or g.survey][resolution]}?" + urlencode({"pos": f"{ra:.5f} {dec:.5f}",
                                                                       "size": f"{size:.1f}"})


def mosaic_tile_centres(g: Giant) -> list[tuple[float, float]]:
    """Centres of the 2 x 2 tiles: offsets of +-0.36 tile north/south and
    east/west of the host, measured on the sky (safe near the pole)."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    c = SkyCoord(g.ra * u.deg, g.dec * u.deg)
    off = 0.36 * MOSAIC_TILE_ARCMIN * u.arcmin
    out = []
    for dn in (-1, 1):
        for de in (-1, 1):
            p = c.directional_offset_by(0 * u.deg, dn * off).directional_offset_by(90 * u.deg, de * off)
            out.append((float(p.ra.deg), float(p.dec.deg)))
    return out


def ps1_url(g: Giant, *, zoom_arcsec: float | None = None) -> str:
    if zoom_arcsec is None:
        size_arcsec = g.size_arcmin * 60
        npx = int(round(min(LEGACY_MAX_PX, size_arcsec / 0.262)))
    else:
        size_arcsec, npx = zoom_arcsec, int(round(zoom_arcsec / 0.262))
    return f"{HIPS2FITS_URL}?" + urlencode({"hips": PS1_HIPS, "width": npx, "height": npx,
                                            "fov": f"{size_arcsec / 3600:.6f}", "projection": "TAN",
                                            "coordsys": "icrs", "ra": f"{g.ra:.6f}", "dec": f"{g.dec:.6f}",
                                            "format": "jpg"})


def nvss_url(g: Giant) -> str:
    pixels = int(round(g.size_arcmin * 60 / 7.5))  # 7.5" pixels: six per NVSS beam
    return f"{SKYVIEW_URL}?" + urlencode({"Survey": "NVSS", "position": f"{g.ra:.5f},{g.dec:.5f}",
                                          "Size": f"{g.size_arcmin / 60:.4f}", "Pixels": pixels,
                                          "Projection": "Sin", "Sampler": "NN", "Return": "FITS"})


def legacy_url(g: Giant, *, zoom_arcsec: float | None = None) -> str:
    if zoom_arcsec is None:
        size_arcsec = g.size_arcmin * 60
        pixscale = max(0.262, size_arcsec / LEGACY_MAX_PX)
    else:
        size_arcsec, pixscale = zoom_arcsec, 0.262
    npx = int(round(size_arcsec / pixscale))
    return f"{LEGACY_URL}?" + urlencode({"ra": f"{g.ra:.5f}", "dec": f"{g.dec:.5f}", "layer": LEGACY_LAYER,
                                         "pixscale": f"{pixscale:.4f}", "width": npx, "height": npx})


#: Close-up of each host: 80" square at the survey's native 0.262"/pixel.
HOST_ZOOM_ARCSEC = 80.0


LOTSS_REF = {"dr2": "Shimwell et al. 2022, A&A 659, A1", "dr3": "Shimwell et al. 2026, arXiv:2602.15949"}


def image_sources() -> list[Source]:
    out = []
    for g in GIANTS:
        rel = g.survey.upper()
        if g.tiles > 1:
            for i, (ra, dec) in enumerate(mosaic_tile_centres(g)):
                out.append(Source(f"images/{g.key}_lotss20_tile{i}.fits",
                                  lotss_url(g, 20, pos=(ra, dec), size_arcmin=MOSAIC_TILE_ARCMIN),
                                  f"LoTSS {rel} 20\" mosaic tile {i + 1}/4 of {g.name}", LOTSS_REF[g.survey],
                                  kind="fits", pause=1.0))
        else:
            for res in (6, 20):
                out.append(Source(f"images/{g.key}_lotss{res}.fits", lotss_url(g, res),
                                  f"LoTSS {rel} {res}\" image of {g.name}", LOTSS_REF[g.survey],
                                  kind="fits", pause=1.0))
        if g.survey == "dr2":
            # the same field in DR3, to check DR2 against its successor
            out.append(Source(f"images/{g.key}_dr3lotss20.fits", lotss_url(g, 20, release="dr3"),
                              f"LoTSS DR3 20\" image of {g.name} (DR2 cross-check)", LOTSS_REF["dr3"],
                              input=False, kind="fits", pause=1.0))
        if g.optical == "ps1":
            out.append(Source(f"images/{g.key}_legacy.jpg", ps1_url(g),
                              f"Pan-STARRS1 colour image of the field of {g.name} (CDS hips2fits)",
                              "Chambers et al. 2016, arXiv:1612.05560", kind="jpeg", pause=1.0))
            out.append(Source(f"images/{g.key}_host.jpg", ps1_url(g, zoom_arcsec=HOST_ZOOM_ARCSEC),
                              f"Pan-STARRS1 close-up of {g.name}'s host galaxy",
                              "Chambers et al. 2016, arXiv:1612.05560", kind="jpeg", pause=1.0))
        else:
            out.append(Source(f"images/{g.key}_legacy.jpg", legacy_url(g),
                              f"Legacy Surveys colour image of the field of {g.name}",
                              "Dey et al. 2019, AJ 157, 168", kind="jpeg", pause=1.0))
            out.append(Source(f"images/{g.key}_host.jpg", legacy_url(g, zoom_arcsec=HOST_ZOOM_ARCSEC),
                              f"Legacy Surveys close-up of {g.name}'s host galaxy",
                              "Dey et al. 2019, AJ 157, 168", kind="jpeg", pause=1.0))
        if "ageing" in g.groups or "gallery" in g.groups:
            out.append(Source(f"images/{g.key}_nvss.fits", nvss_url(g), f"NVSS 1.4 GHz image of {g.name}",
                              "Condon et al. 1998, AJ 115, 1693 (via SkyView)", kind="fits", pause=2.0))
    return out


def all_static_sources() -> list[Source]:
    return CATALOG_SOURCES + vac_sources() + image_sources()
