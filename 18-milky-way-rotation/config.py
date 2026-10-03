"""Project configuration: paths, data sources, and every number we assume.

The analysis measures one thing from the sky, the Doppler velocities of
hydrogen gas. To turn those into a rotation curve and a map, it needs a few
numbers it cannot measure itself: how far the Sun is from the Galactic centre,
how fast the Sun itself moves, and how the survey's velocity frame is defined.
All of them live here, with their sources, so they can be read and changed in
one place.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
HI4PI_DIR = DATA_DIR / "hi4pi"
FIGURE_DIR = PROJECT_DIR / "figures"
PRODUCTS_DIR = DATA_DIR / "products"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, HI4PI_DIR, FIGURE_DIR, PRODUCTS_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# The 21-cm survey: HI4PI
# --------------------------------------------------------------------------

#: HI4PI (HI4PI Collaboration 2016, A&A 594, A116) at the CDS archive. The
#: survey merges the Effelsberg-Bonn HI Survey (north of the celestial
#: equator) with the Parkes Galactic All-Sky Survey (south of it): 16.2'
#: resolution, 1.29 km/s channels, 43 mK noise.
HI4PI_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/594/A116/"
HI4PI_PAPER = "HI4PI Collaboration: Ben Bekhti et al. 2016, A&A 594, A116"

#: The survey is released as 20 x 20 degree cubes in Galactic coordinates
#: (plate-carree projection, 5' pixels, 266 x 266 pixels x 933 channels, 264 MB
#: each). Row "E" is centred on the Galactic plane: 18 cubes, l = 10, 30, ...,
#: 350, each spanning b = -11 to +11 degrees.
HI4PI_ROW_E = [f"CAR_E{i:02d}.fits" for i in range(1, 19)]

#: Seven cubes of row "F" (centred on b = +20, covering b = +9 to +31) at
#: l = 20 to 160. Only the outer-disc map uses them: there the disc warps
#: north, and beyond R ~ 16 kpc the top of the warped layer rises above
#: b = +11 (checked: 20-44% of map cells at azimuths 45-140 deg were cut off
#: with row E alone). The southern warp stays inside row E.
HI4PI_ROW_F_WARP = [f"CAR_F{i:02d}.fits" for i in range(2, 9)]

#: We download only the channels with |v_LSR| <= this, by HTTP byte ranges:
#: each velocity channel is a contiguous 283 kB plane in the file. All Galactic
#: disc gas, including the fast-moving gas near the Galactic centre (up to about
#: +-270 km/s), lies inside it; what is left out is high-velocity clouds and
#: the Magellanic system. It cuts the download from 4.8 GB to 2.8 GB.
V_KEEP_KMS = 350.0

#: 21-cm brightness temperature to column density, optically thin gas:
#: N_HI = C * integral(T_b dv), C in cm^-2 per (K km/s).
C_HI = 1.823e18

#: HI4PI angular resolution (FWHM), degrees.
BEAM_FWHM_DEG = 16.2 / 60


# --------------------------------------------------------------------------
# Galactic constants (the assumed inputs)
# --------------------------------------------------------------------------

#: Reid et al. (2019, ApJ 885, 131), fit A5 to ~200 maser parallaxes and proper
#: motions. R0: Sun to Galactic centre distance; THETA0: circular rotation speed
#: at the Sun; the Sun's own peculiar motion (U toward the Galactic centre, V in
#: the direction of rotation, W toward the north Galactic pole).
REID19_PAPER = "Reid et al. 2019, ApJ 885, 131"
R0_KPC = 8.15
R0_ERR_KPC = 0.15
THETA0_KMS = 236.0
THETA0_ERR_KMS = 7.0
U_SUN_KMS = 10.6
U_SUN_ERR_KMS = 1.2
V_SUN_KMS = 10.7
W_SUN_KMS = 7.6

#: The Sun's full orbital angular velocity, (THETA0 + V_sun) / R0, in km/s/kpc.
#: This is the combination the data pin down best (Reid et al. 2019), because it
#: is the reflex of Sgr A*'s apparent motion across the sky: the black hole sits
#: still at the centre, so the Sun's orbit makes it drift at 6.4 mas/yr.
OMEGA_SUN = 30.32
OMEGA_SUN_ERR = 0.27

#: The *old* IAU standard values (Kerr & Lynden-Bell 1986). Most published
#: terminal-velocity rotation curves, including McClure-Griffiths & Dickey
#: (2007, 2016), are expressed with these, so we also compute ours this way for
#: a like-for-like comparison.
R0_IAU_KPC = 8.5
THETA0_IAU_KMS = 220.0

#: The LSRK velocity frame of the survey is *defined* by the "standard solar
#: motion": 20 km/s toward RA 18h, Dec +30 deg (equinox B1900). It is a
#: convention, not a measurement of the Sun's true motion, and the analysis
#: converts from it explicitly (see galaxy.py).
STANDARD_SOLAR_MOTION_KMS = 20.0
STANDARD_SOLAR_APEX_B1900 = (270.0, 30.0)  # RA, Dec in degrees


# --------------------------------------------------------------------------
# Terminal-velocity measurement
# --------------------------------------------------------------------------

#: Longitude ranges (degrees) for the inner-Galaxy terminal velocities. Inside
#: |l| < 15 deg the Galactic bar drives strongly non-circular orbits; beyond
#: |sin l| ~ 0.95 the terminal velocity shrinks to the size of random motions.
Q1_RANGE = (15.0, 75.0)
Q4_RANGE = (285.0, 345.0)

#: Latitude range averaged before measuring, as in McClure-Griffiths & Dickey.
TERMINAL_B_HALFWIDTH_DEG = 0.5

#: Brightness-temperature thresholds (K). The threshold only finds the edge;
#: the terminal velocity itself comes from fitting the edge's shape, and the
#: result should barely depend on which threshold started the fit.
#: McClure-Griffiths & Dickey used 10 K (2007) and 20 K (2016).
TERMINAL_THRESHOLDS_K = (10.0, 20.0)
TERMINAL_THRESHOLD_K = 10.0

#: The edge fit uses velocities from this far inside the threshold crossing
#: (km/s) ... to this far beyond it. 25 km/s is the narrowest window inside
#: which every fitted edge centre lands at every longitude in both quadrants;
#: narrower windows clip bright edges, whose half-power point sits deep
#: inside the 10 K crossing. The result's sensitivity to it is reported.
TERMINAL_W_IN_KMS = 25.0
TERMINAL_W_OUT_KMS = 45.0
TERMINAL_W_IN_VARIANTS = (20.0, 25.0, 30.0)

#: A longitude's terminal velocity counts as reliable only if the 12 method
#: variants (2 thresholds x 3 windows x mean/median over latitude) agree to
#: within this rms (km/s). Where they do not, the edge itself is ambiguous:
#: several emission components, absorption against a bright continuum source,
#: or the bar's non-circular gas near l = 343.
TERMINAL_MAX_METHOD_SPREAD = 3.0


# --------------------------------------------------------------------------
# Reference tables (small; fetched alongside the survey)
# --------------------------------------------------------------------------

SPARC_URL = "https://astroweb.cwru.edu/SPARC/"


@dataclass(frozen=True)
class Source:
    """One small file we download, and why."""

    filename: str
    url: str
    what: str
    reference: str
    #: ``True`` if the analysis is built on it; ``False`` if it is a published
    #: result used only to grade our own at the end.
    input: bool = True


SOURCES: list[Source] = [
    Source(
        "mgd16_q1_terminal.dat",
        "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/831/124/table1.dat",
        "Published first-quadrant HI terminal velocities, 748 longitudes (VLA Galactic "
        "Plane Survey, 1' resolution)",
        "McClure-Griffiths & Dickey 2016, ApJ 831, 124", input=False,
    ),
    Source(
        "reid19_masers.dat",
        "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/885/131/table1.dat",
        "Parallaxes, proper motions and velocities of 199 masers in young massive-star "
        "regions, with spiral-arm assignments",
        REID19_PAPER, input=False,
    ),
    Source(
        "SPARC_Lelli2016c.mrt", SPARC_URL + "SPARC_Lelli2016c.mrt",
        "SPARC Table 1: 175 disc galaxies (distances, luminosities, flat rotation speeds)",
        "Lelli, McGaugh & Schombert 2016, AJ 152, 157", input=False,
    ),
    Source(
        "MassModels_Lelli2016c.mrt", SPARC_URL + "MassModels_Lelli2016c.mrt",
        "SPARC Table 2: rotation curves of the 175 galaxies",
        "Lelli, McGaugh & Schombert 2016, AJ 152, 157", input=False,
    ),
]
