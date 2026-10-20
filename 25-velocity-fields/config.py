"""Project configuration: paths, data sources, the galaxies, and every assumed number.

Everything the analysis *assumes* rather than measures lives here, so the
assumptions can be read in one place. Published results (de Blok et al. 2008
and friends) also live here or in files fetched alongside the data, but the
analysis never starts from them: they are read only at the end, to grade our
own numbers.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
THINGS_DIR = DATA_DIR / "things"
REF_DIR = DATA_DIR / "references"
MANGA_DIR = DATA_DIR / "manga"
MANGA_BINS_DIR = MANGA_DIR / "bins"
MANGA_MAPS_DIR = MANGA_DIR / "maps"
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, THINGS_DIR, REF_DIR, MANGA_DIR, MANGA_BINS_DIR, MANGA_MAPS_DIR,
           DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-project-atlas/J1 (hobby research; python-urllib)"

#: Speed of light, km/s (exact).
C_KMS = 299_792.458
#: Rest frequency of the 21-cm hyperfine line of neutral hydrogen, Hz.
HI_REST_HZ = 1.420405752e9
#: Newton's constant in kpc (km/s)^2 / Msun.
G_KPC = 4.300917e-6


# --------------------------------------------------------------------------
# THINGS: The HI Nearby Galaxy Survey (Walter et al. 2008)
# --------------------------------------------------------------------------

#: The THINGS archive at MPIA. The old address (www.mpia.de/THINGS) redirects here.
THINGS_URL = "https://things.www3.mpia.de/Data_files/"
THINGS_PAPER = "Walter et al. 2008, AJ 136, 2563"


@dataclass(frozen=True)
class Galaxy:
    """A THINGS galaxy we turn from a data cube into a rotation curve."""

    #: Display name, e.g. "NGC 3198".
    name: str
    #: File-name stem in the THINGS archive, e.g. "NGC_3198".
    things: str
    #: Name in SPARC (project D4), e.g. "NGC3198".
    sparc: str
    #: What it is, and why it is here.
    why: str
    #: Inside this radius (arcsec) the position angle and inclination are not
    #: fitted but carried inward from the first ring outside it. The default
    #: (None) is two beams, inside which the geometry is not resolved; a bar
    #: needs more, because gas on bar orbits mimics a tilted disc.
    geometry_fixed_inside_arcsec: float | None = None
    #: Why that radius, when it is not the default.
    geometry_note: str = ""

    @property
    def key(self) -> str:
        return self.name.replace(" ", "")

    def cube_file(self, weighting: str = "NA") -> str:
        return f"{self.things}_{weighting}_CUBE_THINGS.FITS"

    def moment_file(self, moment: int, weighting: str = "NA") -> str:
        return f"{self.things}_{weighting}_MOM{moment}_THINGS.FITS"


GALAXIES: list[Galaxy] = [
    Galaxy("NGC 3198", "NGC_3198", "NGC3198",
           "The classic flat rotation curve: van Albada et al. (1985) used it to show "
           "that no disc of stars can do without a dark halo. Our flagship."),
    Galaxy("NGC 2403", "NGC_2403", "NGC2403",
           "A nearby, regular spiral imaged at 1-arcsecond pixels; its gas streams "
           "along the spiral arms."),
    Galaxy("NGC 2903", "NGC_2903", "NGC2903",
           "A barred spiral. The bar drives gas on elongated orbits, so the inner "
           "velocity field is not circular rotation at all.",
           geometry_fixed_inside_arcsec=110.0,
           geometry_note="Free fits find the inclination falling 13 degrees between 77 and "
           "119 arcsec, the outer edge of the bright spiral arms that wrap the bar, and "
           "steady beyond; inside it the bar's streaming motions masquerade as a changing "
           "tilt. de Blok et al. (2008) made the same choice."),
    Galaxy("DDO 154", "DDO154", "DDO154",
           "A gas-rich dwarf whose 21-cm disc extends far beyond its stars: the test "
           "case for the cusp-core problem."),
]

FLAGSHIP = "NGC 3198"


def galaxy_by_name(name: str) -> Galaxy:
    key = name.replace(" ", "").replace("_", "").lower()
    for g in GALAXIES:
        if g.key.lower() == key:
            return g
    raise KeyError(f"{name!r} is not one of {[g.name for g in GALAXIES]}")


# --------------------------------------------------------------------------
# Making velocity fields from the cubes
# --------------------------------------------------------------------------

#: The quality criteria of de Blok et al. (2008, Sect. 3.4) for a fitted
#: profile to enter the velocity field: its peak must exceed this many times
#: the channel noise, and its (equivalent Gaussian) dispersion must exceed the
#: channel spacing. A third criterion masks positions whose integrated HI is
#: below this many times the noise of the integrated map.
PEAK_SNR_MIN = 3.0
MOM0_SNR_MIN = 3.0

#: Channels at each end of a cube treated as line-free when no emission is
#: found there (checked per cube, see things.py).
NOISE_EDGE_CHANNELS = 6


# --------------------------------------------------------------------------
# Tilted rings
# --------------------------------------------------------------------------

#: Rings are spaced at half the beam (two points per beam), as de Blok et al.
#: (2008) did; consecutive rings therefore share information.
RINGS_PER_BEAM = 2

#: Weight each pixel by |cos(theta)|^WEIGHT_POWER, theta the azimuth in the
#: disc. Pixels near the minor axis carry almost no rotation signal. de Blok et
#: al. (2008) used |cos theta|, i.e. power 1.
WEIGHT_POWER = 1

#: Box-car width (in rings) used to smooth the free PA and inclination profiles
#: before the final rotation-curve fit. de Blok et al. used 3 or 5.
GEOMETRY_SMOOTH_RINGS = 5

#: Minimum number of velocity-field pixels a ring needs for a fit, as a
#: multiple of the number of beams it covers. Rings at the ragged outer edge
#: with only a few clumps of gas are dropped rather than over-interpreted.
MIN_BEAMS_PER_RING = 4.0

#: Highest order of the harmonic decomposition (Schoenmakers 1999;
#: Trachternach et al. 2008 found 3 enough).
HARMONIC_ORDER = 3


# --------------------------------------------------------------------------
# Mass models
# --------------------------------------------------------------------------

#: The mass in gas per unit HI mass: helium and metals ride along with the
#: hydrogen that the 21-cm line sees. 1.4 is the value Oh et al. (2011) and
#: de Blok et al. (2008) used; SPARC (D4) used 1.33.
GAS_PER_HI = 1.4

#: HI column density per unit 21-cm brightness: N_HI = 1.823e18 cm^-2 per
#: (K km/s), for optically thin gas.
NHI_PER_K_KMS = 1.823e18

#: Vertical structure assumed for the beam-smearing model galaxies: an HI layer
#: of this Gaussian scale height (kpc) and velocity dispersion (km/s). Values
#: typical of THINGS discs (de Blok et al. 2008 used 0.1 kpc and 8 km/s); the
#: test is insensitive to them.
MODEL_SCALE_HEIGHT_KPC = 0.1
MODEL_DISPERSION_KMS = 8.0


# --------------------------------------------------------------------------
# MaNGA (SDSS-IV DR17)
# --------------------------------------------------------------------------

MANGA_DAP_URL = "https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/"
MANGA_DRP_URL = "https://data.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/"
MANGA_DAPALL = "dapall-v3_1_1-3.1.0.fits"
MANGA_DRPALL = "drpall-v3_1_1.fits"
#: Voronoi-binned to S/N 10 for both stars and gas: one velocity per bin for
#: each, measured on the same bins, which is what a stars-vs-gas comparison
#: wants.
MANGA_DAPTYPE = "VOR10-MILESHC-MASTARSSP"
MANGA_PAPER = "Bundy et al. 2015, ApJ 798, 7; Westfall et al. 2019, AJ 158, 231 (DAP)"

#: Parallel connections to the SDSS archive for the survey-wide download.
MANGA_WORKERS = 8

#: The kinematic position angle is fitted to bins within this many effective
#: radii, whose spectra reach this median S/N per pixel -- the choices of Zhou et
#: al. (2022), so our misaligned fractions can be compared with theirs.
PA_FIT_MAX_RE = 1.5
PA_FIT_MIN_SNR = 3.0
#: A position angle counts as measured if its 1-sigma error is at most this.
PA_MAX_ERROR_DEG = 60.0
#: Misaligned: stars and gas rotate about axes more than 30 degrees apart.
#: Counter-rotating: more than 150 degrees apart (anti-parallel within 30).
MISALIGNED_DEG = 30.0
COUNTER_ROTATING_DEG = 150.0

#: For the asymmetric-drift showcase: the intrinsic thickness assumed when
#: turning a disc's apparent axis ratio into an inclination (Hubble 1926).
INTRINSIC_AXIS_RATIO = 0.2
#: The ratio of vertical to radial velocity dispersion of disc stars, used to
#: turn MaNGA's line-of-sight dispersion into the radial one. 0.6 +- 0.15 in
#: the solar neighbourhood and in other discs (e.g. Gerssen & Shapiro Griffin
#: 2012; Martinsson et al. 2013).
SIGMA_Z_OVER_R = 0.6
