"""Project configuration: paths, the Galaxy model, every selection number, and
the published values we grade ourselves against.

H3 has three parts, each a different way of reading the Milky Way's history out
of Gaia's six-dimensional catalogue (position on the sky, distance, two proper
motions, and radial velocity):

1. **GD-1**, a thin stream of stars torn out of a globular cluster that no
   longer exists. We are told *where* on the sky it lies -- it was discovered
   in 2006 -- but not how it moves or how far away it is. Both are found in
   the data.
2. **Gaia-Sausage-Enceladus**, the remains of a galaxy the Milky Way swallowed
   ~10 Gyr ago, found in the velocities and orbital energies of nearby halo
   stars -- once with a velocity cut, as the discovery papers did, and once
   with chemistry alone, which avoids the circular argument that a velocity cut
   invites.
3. **The phase spiral**, a snail-shell pattern in the vertical positions and
   velocities of stars near the Sun: the disc is still ringing from a
   disturbance a few hundred million years ago. We date it.

What the configuration deliberately does NOT contain is GD-1's proper motion or
distance, or any answer the analysis is supposed to find. Published values live
in the ``REFERENCE`` block at the bottom and are used only to grade results.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
GAIA_DIR = DATA_DIR / "gaia"
CHUNK_DIR = GAIA_DIR / "chunks"
AUX_DIR = DATA_DIR / "aux"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, GAIA_DIR, CHUNK_DIR, AUX_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# The Sun's place in the Galaxy
# --------------------------------------------------------------------------
# astropy's v4.0 Galactocentric defaults, so every number here can be checked
# against a standard library in one line (test_physics.py does).

#: Distance from the Sun to the Galactic Centre (GRAVITY Collaboration 2018,
#: from the orbit of the star S2 -- project H1 measures this itself).
R0_KPC = 8.122
#: Height of the Sun above the Galactic midplane (Bennett & Bovy 2019).
Z_SUN_KPC = 0.0208
#: The Sun's velocity in the Galactocentric frame, km/s: towards the Galactic
#: Centre, in the direction of rotation, and towards the North Galactic Pole
#: (Drimmel & Poggio 2018: Sgr A*'s proper motion, Reid & Brunthaler 2004,
#: combined with R0).
V_SUN_KMS = (12.9, 245.6, 7.78)
#: The Sun's motion relative to the Local Standard of Rest -- the velocity of
#: a circular orbit at the Sun's radius (Schoenrich, Binney & Dehnen 2010).
#: Used only to measure how fast a star moves relative to the rotating disc,
#: the quantity the halo selection of Helmi et al. (2018) cuts on.
V_PECULIAR_KMS = (11.1, 12.24, 7.25)

#: 1 mas/yr at 1 kpc, in km/s.
K_KMS = 4.740470463533348

# --------------------------------------------------------------------------
# Gaia DR3 access
# --------------------------------------------------------------------------

GAIA_SYNC_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
#: Fallback mirror of Gaia DR3 at CDS (used by projects B2/B3 when ESA's
#: archive was down).
VIZIER_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"

#: Every Gaia DR3 source_id encodes its sky position: source_id // 2**35 is
#: the star's HEALPix pixel at level 12, nested ordering. So a *range* of
#: source_ids is a patch of sky, and the table is stored in that order -- a
#: range query is the fastest way to cut the sky into pieces. Level 1 splits
#: the sky into 48 patches of ~860 square degrees each.
HEALPIX_BASE_LEVEL = 1
#: A patch whose query fails (ESA's synchronous endpoint gives up after about
#: a minute) is split into its four children and retried, down to this level.
HEALPIX_MAX_LEVEL = 4

#: A random 2% of all 6D stars, drawn with Gaia's own random_index column (a
#: random permutation of every source in DR3), for context: the disc that the
#: halo selections are cut out of.
RANDOM_FRACTION = 0.02
N_DR3_SOURCES = 1_811_709_771

# --------------------------------------------------------------------------
# The 6D sample shared by parts 2 and 3
# --------------------------------------------------------------------------

#: Stars with a Gaia radial velocity and a parallax measured to better than
#: 20%. Distance is then simply 1/parallax -- the choice made by Antoja et al.
#: (2018) and Helmi et al. (2018), whose results we check against. Gaia DR3
#: has 29.9 million such stars (counted on the archive, 30 Sep 2026).
PARALLAX_OVER_ERROR_6D = 5.0

#: Quality: renormalised unit weight error. Above ~1.4 the single-star
#: astrometric model fits badly, usually an unresolved binary (Lindegren et al.
#: 2021). Applied locally, so its effect can be shown.
RUWE_MAX = 1.4

# --- Part 3: the phase spiral ---------------------------------------------

#: Stars within this distance of the Sun's Galactocentric radius. Antoja et
#: al. (2018) used 8.24 < R < 8.44 kpc for their R0 = 8.34 kpc: a ring
#: 0.2 kpc wide centred on the Sun. Same width here, centred on our R0.
RING_HALF_WIDTH_KPC = 0.1
#: The ring is downloaded a little wider than it is used, so that a
#: different R0 can be tried without re-fetching.
RING_FETCH_HALF_WIDTH_KPC = 0.15

# --- Part 2: the stellar halo ----------------------------------------------

#: Stars moving at least this fast relative to the Local Standard of Rest are
#: downloaded in full. The halo definition we test (Helmi et al. 2018) is
#: 210 km/s; fetching from 150 keeps the fast tail of the thick disc too, so
#: the boundary itself can be moved and its effect seen.
HALO_FETCH_SPEED_KMS = 150.0
#: The halo definition of Helmi et al. (2018): |v - v_LSR| > 210 km/s.
HALO_SPEED_KMS = 210.0
#: Helmi et al. (2018) Fig. 1 uses stars within 2.5 kpc (parallax > 0.4 mas)
#: "to limit the impact of velocity gradients".
HALO_MAX_DISTANCE_KPC = 2.5

#: Gaia's own chemistry: the GSP-Spec pipeline's metallicities from RVS
#: spectra (Recio-Blanco et al. 2023). Stars more metal-poor than this (raw,
#: uncalibrated [M/H]) are downloaded with no cut on velocity at all -- the
#: non-circular route to the halo that Belokurov et al. (2018) took.
GSPSPEC_FETCH_MH_MAX = -0.7

#: The first 13 characters of flags_gspspec describe how trustworthy the
#: atmospheric parameters are: 1-3 line broadening, 4-6 radial-velocity
#: correction, 7 flux noise, 8 extrapolation beyond the model grid, 9-12 bad
#: or missing flux and emission, 13 cool K/M giants. 0 = best, 9 = worst.
#: The "medium-quality" sample of Recio-Blanco et al. (2023, Sect. 9.1):
#: every flag <= 1, except the flux-noise flag (relaxed to <= 3) and the cool
#: giant flag (which must be 0). Maximum allowed value per flag:
GSPSPEC_MEDIUM_MAX = (1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 0)

#: GSP-Spec calibration polynomials in uncalibrated log(g) (Recio-Blanco et
#: al. 2023, Tables 3 and 4). [M/H] is corrected as [M/H] + sum p_i logg^i.
GSPSPEC_MH_CALIB = (0.274, -0.1373, -0.0050, 0.0048)
#: [alpha/Fe] calibration, valid for 1.01 < log(g) < 4.85 (extrapol flag 0).
GSPSPEC_ALPHA_CALIB = (-0.5809, 0.7018, -0.2402, 0.0239)
GSPSPEC_ALPHA_LOGG_RANGE = (1.01, 4.85)

#: Metallicity bins for the anisotropy test. Belokurov et al. (2018) used
#: -2.30, -2.00, -1.66, -1.33, -1.00 as bin edges.
FEH_BINS = (-3.0, -2.0, -1.66, -1.33, -1.0, -0.7, -0.4, 0.0, 0.5)

# --------------------------------------------------------------------------
# The Milky Way's gravity: two published potentials
# --------------------------------------------------------------------------
# An orbit's energy depends on the assumed potential, so published selection
# boxes in energy only mean something in the potential they were drawn in.
# We compute energies in both potentials the literature boxes use.


@dataclass(frozen=True)
class Helmi17Potential:
    """The potential of Helmi et al. (2017), used by Helmi et al. (2018) and
    Koppelman et al. (2019): Miyamoto-Nagai disc + Hernquist bulge + NFW halo,
    parameters as listed in Koppelman et al. (2019), Sect. 3."""

    disc_mass_msun: float = 9.3e10
    disc_a_kpc: float = 6.5
    disc_b_kpc: float = 0.26
    bulge_mass_msun: float = 3.0e10
    bulge_c_kpc: float = 0.7
    halo_mass_msun: float = 1.0e12  # virial mass
    halo_rs_kpc: float = 21.5
    halo_concentration: float = 12.0


HELMI17 = Helmi17Potential()

#: Newton's G in kpc (km/s)^2 / Msun (IAU 2015 nominal solar mass parameter).
G_KPC_KMS2_MSUN = 4.300917270038e-06

# --------------------------------------------------------------------------
# Part 1: GD-1
# --------------------------------------------------------------------------

#: The stream frame of Koposov et al. (2010), Appendix: phi1 runs along the
#: stream, phi2 across it. This is the stream's *position* on the sky, known
#: since its discovery (Grillmair & Dionatos 2006) -- the one thing about
#: GD-1 the pipeline is told.
GD1_ROTATION = (
    (-0.4776303088, -0.1738432154, 0.8611897727),
    (0.510844589, -0.8524449229, 0.111245042),
    (0.7147776536, 0.4930681392, 0.4959603976),
)

#: The field downloaded. At phi1 = -90 deg the stream runs into the Galactic
#: plane (b = 13 deg, and 4 deg at -100); the phi1 < -80 strip alone held
#: 1.06 million stars, more than a third of the whole field, and PWB2018 see
#: the stream only at b > 20 deg. So the field stops at -80.
GD1_PHI1_RANGE = (-80.0, 15.0)
GD1_PHI2_RANGE = (-8.0, 6.0)
#: The field is fetched in pieces this long in phi1, each through a circular
#: sky search (fast on ESA's archive, whereas polygon searches time out).
GD1_QUERY_SEGMENT_DEG = 5.0

#: Archive-side cuts. GD-1 is ~8-10 kpc away (parallax ~0.1 mas), so
#: parallax < 1 mas is the cut Price-Whelan & Bonaca (2018) used to drop
#: nearby stars. BP-RP < 1.6 drops red dwarfs -- the bulk of the faint
#: foreground -- while keeping GD-1's whole main sequence, turn-off and giant
#: branch (dereddened BP-RP 0.5-1.3 for an old, metal-poor population).
GD1_G_MAX = 20.5
GD1_PARALLAX_MAX_MAS = 1.0
GD1_BP_RP_MAX = 1.6

#: Bin size of the unfiltered "before" star-count map, degrees.
GD1_COUNT_MAP_BIN_DEG = 0.25

#: On-stream and off-stream bands in phi2, degrees, for the blind search for
#: the stream's proper motion. The off-stream control skips the regions of
#: the spur and blob reported by PWB2018.
GD1_ON_STREAM_PHI2 = 0.75
GD1_OFF_STREAM_PHI2 = (2.0, 5.5)
#: Pieces of the stream searched separately, degrees in phi1.
GD1_SEARCH_SEGMENT_DEG = 10.0
#: Proper-motion histogram bin and search window for the knot, mas/yr.
GD1_PM_BIN = 0.25
GD1_PM_WINDOW = ((-20.0, 5.0), (-10.0, 10.0))

# --- M92: an old, metal-poor template ------------------------------------

@dataclass(frozen=True)
class TemplateCluster:
    """A globular cluster whose colour-magnitude diagram serves as the filter
    for GD-1's stars: the same kind of population (old, metal-poor), measured
    with the same instrument, so no stellar-model isochrone is needed."""

    name: str
    ra_deg: float
    dec_deg: float
    #: Cone radius fetched, degrees. M92's tidal radius is 15.2 arcmin (Harris
    #: 1996, 2010 edition).
    radius_deg: float
    #: The crowded core is skipped: Gaia's BP and RP photometry blends there.
    core_exclude_deg: float
    distance_kpc: float
    distance_err_kpc: float
    distance_source: str
    feh: float
    feh_source: str
    ebv: float


M92 = TemplateCluster(
    name="M92 (NGC 6341)",
    ra_deg=259.28079,
    dec_deg=43.13594,
    radius_deg=0.30,
    core_exclude_deg=0.04,
    distance_kpc=8.501,
    distance_err_kpc=0.071,
    distance_source="Baumgardt & Vasiliev (2021), MNRAS 505, 5957",
    feh=-2.31,
    feh_source="Harris (1996, 2010 edition)",
    ebv=0.02,
)

# --- Dust ------------------------------------------------------------------

#: Schlegel, Finkbeiner & Davis (1998) E(B-V) map, HEALPix version from NASA's
#: LAMBDA archive (25 MB).
SFD_URL = "https://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits"
#: Schlafly & Finkbeiner (2011): SFD overestimates reddening by 14%.
SFD_RECALIBRATION = 0.86
#: A_V / E(B-V) for the Milky Way's average dust.
R_V = 3.1

# --------------------------------------------------------------------------
# Reference values -- used only to grade the results, never as inputs
# --------------------------------------------------------------------------

REFERENCE = {
    "antoja2018": {
        "citation": "Antoja et al. (2018), Nature 561, 360",
        "n_ring_stars_dr2": 935_590,
        "n_6d_dr2": 6_376_803,
        "ring_kpc": (8.24, 8.44),
        "r0_kpc": 8.34,
        "perturbation_myr": 500.0,
        "perturbation_range_myr": (300.0, 900.0),
    },
    "pwb2018": {
        "citation": "Price-Whelan & Bonaca (2018), ApJL 863, L20",
        # Reflex-corrected proper-motion box that selects GD-1 (their Fig. 1).
        "pm_box_reflex": ((-9.0, -4.5), (-1.7, 1.0)),
        # Mean proper motion at phi1 in (-18, -10) deg, not reflex corrected.
        "pm_radec_at_m14": (-7.78, -7.85),
        "radec_at_m14": (177.01, 53.99),
        "spur_phi1": (-36.0, -30.0),
        "spur_phi2": (0.9, 1.7),
        "blob_phi1": (-18.0, -12.0),
        "control_phi1": (-51.0, -45.0),
        "progenitor_candidate": (-13.5, -0.5),
        "gap_phi1": -20.0,
        "stream_half_width_deg": 0.75,
        "extent_phi1": (-90.0, 10.0),
    },
    "bonaca2019": {
        "citation": "Bonaca et al. (2019), ApJ 880, 38",
        # Two significant gaps, G-20 and G-40; the spur extends ~10 deg from
        # G-40. Their model gap profile is centred on -40.5 deg, 8.5 deg wide.
        "gaps_phi1": (-40.0, -20.0),
        "g40_centre": -40.5,
        "g40_width": 8.5,
        "perturber_mass_msun": (1e6, 1e8),
    },
    "koposov2010": {
        "citation": "Koposov, Rix & Hogg (2010), ApJ 712, 260",
        # Table 3: heliocentric distance along the stream.
        "distance_phi1": (-55.0, -45.0, -35.0, -25.0, -15.0, 0.0),
        "distance_kpc": (7.20, 7.59, 7.83, 8.69, 8.91, 9.86),
        "distance_err_kpc": (0.30, 0.40, 0.30, 0.40, 0.40, 0.50),
        # Table 2: stream track phi2(phi1).
        "track_phi1": (-60.0, -56.0, -54.0, -48.0, -44.0, -40.0, -36.0, -34.0, -32.0,
                       -30.0, -28.0, -24.0, -22.0, -18.0, -12.0, -2.0),
        "track_phi2": (-0.64, -0.89, -0.45, -0.08, 0.01, -0.00, 0.04, 0.06, 0.04,
                       0.08, 0.03, 0.06, 0.06, -0.05, -0.29, -0.87),
        "track_phi2_err": (0.15, 0.27, 0.15, 0.13, 0.14, 0.09, 0.10, 0.13, 0.06,
                           0.10, 0.12, 0.05, 0.13, 0.11, 0.16, 0.07),
        # Table 1: radial velocities of spectroscopic members (phi1, v_r, err).
        "rv_phi1": (-45.23, -43.17, -39.54, -39.25, -37.95, -37.96, -35.49, -35.27,
                    -34.92, -34.74, -33.74, -32.90, -32.25, -29.95, -26.61, -25.45,
                    -24.86, -21.21, -14.47, -13.73, -13.02, -12.68, -12.55),
        "rv_kms": (28.8, 29.3, 2.9, -5.2, 1.1, -11.7, -50.4, -30.9, -35.3, -30.9,
                   -74.3, -71.5, -71.5, -92.7, -114.2, -67.8, -111.2, -144.4,
                   -179.0, -191.4, -162.9, -217.2, -172.2),
        "rv_err_kms": (6.9, 10.2, 8.7, 6.5, 5.6, 11.2, 5.2, 12.8, 7.5, 9.2, 9.8,
                       9.6, 9.2, 8.7, 7.3, 7.1, 17.8, 10.5, 10.0, 7.5, 9.6, 10.7,
                       6.6),
    },
    "helmi2018": {
        "citation": "Helmi et al. (2018), Nature 563, 85",
        "halo_speed_kms": 210.0,
        "max_distance_kpc": 2.5,
        # Gaia-Enceladus in the Helmi17 potential.
        "gse_lz": (-1500.0, 150.0),
        "gse_energy_min": -1.8e5,
        "gse_feh_peak": -1.6,
    },
    "myeong2019": {
        "citation": "Myeong et al. (2019), MNRAS 488, 1235",
        # Metallicity distribution peaks: the Sausage (GSE) and Sequoia.
        "gse_feh_peak": -1.3,
        "sequoia_feh_peak": -1.6,
    },
    "feuillet2021": {
        "citation": "Feuillet et al. (2021), MNRAS 508, 1489",
        # GSE selected by radial action and L_z: a narrow metallicity
        # distribution centred at [Fe/H] = -1.17, dispersion 0.34 dex, which
        # they note is 0.1-0.3 dex more metal-rich than earlier estimates.
        "gse_feh_mean": -1.17,
        "gse_feh_dispersion": 0.34,
    },
    "belokurov2018": {
        "citation": "Belokurov et al. (2018), MNRAS 478, 611",
        # Orbital anisotropy of the halo: beta ~ 0.9 for [Fe/H] > -1.7,
        # 0.2 < beta < 0.4 at the low-metallicity end.
        "beta_metal_rich_halo": 0.9,
        "beta_metal_poor_halo": (0.2, 0.4),
        "feh_transition": -1.7,
    },
    "koppelman2019": {
        "citation": "Koppelman et al. (2019), A&A 625, A5",
        # Helmi streams in angular momentum (box A; box B is 1600-3200,
        # 750-1700).
        "hs_lperp": (1750.0, 2600.0),
        "hs_lz": (1000.0, 1500.0),
    },
    "koppelman2019b": {
        "citation": "Koppelman et al. (2019b), A&A 631, L9",
        # In the McMillan (2017) potential; energy in 1e5 km^2/s^2, and
        # circularity = Lz / Lz(circular orbit of the same energy).
        "gse": {"energy": (-1.5, -1.1), "circ": (-0.2, 0.2)},
        "sequoia": {"energy": (-1.35, -1.0), "circ": (-0.65, -0.4)},
        "thamnos1": {"energy": (-1.65, -1.45), "circ": (-1.0, -0.75)},
        "thamnos2": {"energy": (-1.8, -1.6), "circ": (-0.75, -0.4)},
        "helmi_streams": {"lperp": (1600.0, 3200.0), "lz": (1000.0, 1500.0)},
    },
}
