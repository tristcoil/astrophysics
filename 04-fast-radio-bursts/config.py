"""Project configuration: paths, pinned data sources, and every analysis choice.

Every number that decides which bursts enter which plot lives here, with the
reason it has the value it does. Changing a cut means changing it in exactly one
place; the figures and the numbers in the README follow.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
DYNSPEC_DIR = DATA_DIR / "dynamic_spectra"
FRB_REPO_DIR = DATA_DIR / "frbs_repo"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR, DYNSPEC_DIR, FRB_REPO_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: The Second CHIME/FRB Catalog (CHIME/FRB Collaboration 2026, arXiv:2601.09399):
#: 4539 bursts from 3641 sources seen between 2018 July 25 and 2023 September 15.
#: It supersedes Catalog 1 (536 bursts) -- 533 of those were reprocessed into it.
#: Hosted at the Canadian Advanced Network for Astronomical Research (CANFAR).
CHIME_DOI = "10.11570/25.0066"
_CANFAR = (
    "https://www.canfar.net/storage/vault/file/AstroDataCitationDOI/"
    "CISTI.CANFAR/25.0066/data"
)
CHIME_CATALOG_URL = f"{_CANFAR}/table/chimefrbcat2.csv"
CHIME_CATALOG_PATH = DATA_DIR / "chimefrbcat2.csv"
#: Per-burst Stokes-I dynamic spectra (~10-20 MB each), 16384 channels x ~1 ms.
CHIME_DYNSPEC_URL = f"{_CANFAR}/dynamic_spectra/data/{{tns}}_stokesi_dynamic_spectrum.h5"

#: The FRB community's open compilation of localized bursts and their host
#: galaxies (github.com/FRBs/FRB, maintained by Prochaska et al.). Pinned to one
#: commit so the host sample cannot change underneath the results.
FRB_REPO_COMMIT = "996fcda9b0b22431e3171208b4c8e1cf2798e823"  # 2026-05-06
_FRB_RAW = f"https://raw.githubusercontent.com/FRBs/FRB/{FRB_REPO_COMMIT}/frb/data"
FRB_REPO_FILES = {
    # Localized FRBs: position, DM, redshift, repeater flag, references.
    "FRBs_base.csv": f"{_FRB_RAW}/FRBs/FRBs_base.csv",
    # Host galaxies: PATH association probability P(O|x) and redshift.
    "public_hosts.csv": f"{_FRB_RAW}/Galaxies/public_hosts.csv",
    # The NE2001 Galactic electron model integrated through the whole Galaxy,
    # as a HEALPix map (NSIDE 128, Galactic coordinates, RING ordering). The
    # file's header does not say which frame it is in; analysis.py checks it
    # against the CHIME catalogue's own NE2001 values before using it.
    "ne2001_dm_healpix_map_128.fits": f"{_FRB_RAW}/DM/ne2001_dm_healpix_map_128.fits",
}

#: The ATNF Pulsar Catalogue (Manchester et al. 2005, AJ 129, 1993), used as a
#: control sample of radio sources known to lie *inside* the Milky Way (plus a
#: few dozen in the Magellanic Clouds). Queried through its web form for just
#: the columns we need; the catalogue version is read from the response.
ATNF_URL = "https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php"
ATNF_PARAMS = {
    "JName": "JName", "RaJD": "RaJD", "DecJD": "DecJD", "Gl": "Gl", "Gb": "Gb",
    "DM": "DM", "Assoc": "Assoc", "Type": "Type",
    "startUserDefined": "true", "sort_attr": "jname", "sort_order": "asc",
    "condition": "", "pulsar_names": "", "ephemeris": "short",
    "coords_unit": "raj/decj", "radius": "", "coords_1": "", "coords_2": "",
    "style": "Short without errors", "no_value": "*", "fsize": "3",
    "x_axis": "", "x_scale": "linear", "y_axis": "", "y_scale": "linear",
    "state": "query", "table_bottom.x": "50", "table_bottom.y": "18",
}
ATNF_PATH = DATA_DIR / "atnf_pulsars.csv"


# ---------------------------------------------------------------------------
# Example bursts for the dynamic-spectrum figures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ExampleBurst:
    """One burst whose full dynamic spectrum we download and show."""

    tns: str
    role: str
    why: str


#: Chosen for clarity, not at random, so they illustrate rather than prove. The
#: selection rule for each is stated; the statistics in the figures that sit next
#: to them come from the whole catalogue, not from these.
EXAMPLE_BURSTS: list[ExampleBurst] = [
    ExampleBurst(
        "FRB20211005A", "one-off",
        "Highest-S/N burst in Catalog 2 that is single-component, unscattered, "
        "and fills the band; never seen to repeat. DM 226.",
    ),
    ExampleBurst(
        "FRB20200711F", "one-off",
        "Third-highest by the same rule, picked over the second for its higher "
        "DM (528), which doubles the dispersion sweep next to FRB20211005A.",
    ),
    ExampleBurst(
        "FRB20181222A", "repeater",
        "A bright (S/N 230) burst from FRB 20180916B, the repeater with a 16.3-day "
        "activity cycle, with five sub-bursts drifting downward in frequency.",
    ),
    ExampleBurst(
        "FRB20211114A", "repeater",
        "A burst from repeater FRB 20190117A: four sub-bursts drifting downward "
        "in frequency, the 'sad trombone' pattern seen in many repeaters.",
    ),
    ExampleBurst(
        "FRB20230905D", "scattered",
        "Highest-S/N one-off burst in Catalog 2 with scattering above 5 ms "
        "(84 ms at 400 MHz): the clearest view of multipath smearing.",
    ),
]

#: The two one-offs whose dispersion sweeps are drawn together in figure 02.
DISPERSION_BURSTS = ("FRB20211005A", "FRB20200711F")
#: Burst used to show what a wrong DM does.
WRONG_DM_BURST = "FRB20211005A"
#: Fractional DM error for the "wrong DM" panels. At DM 226, 1% leaves 44 ms of
#: residual sweep across the band -- inside the public data window (+/-77 ms).
WRONG_DM_FRACTION = 0.01
#: Downsampling of the 16384 x 24.4 kHz channels for display: 64 -> 1.56 MHz.
DISPLAY_FREQ_DOWNSAMPLE = 64


# ---------------------------------------------------------------------------
# CHIME catalogue cuts
# ---------------------------------------------------------------------------

#: Flags whose events CHIME says to exclude from population studies: operation
#: during commissioning or low-sensitivity days, far-sidelobe detections, and
#: low-S/N events recovered by citizen scientists. Applied to the first burst of
#: each source, following the catalogue paper's convention.
POPULATION_FLAGS = ("excluded_flag", "sidelobe_flag", "citizen_science_flag")
#: For the morphology comparison we additionally drop fits the catalogue marks
#: as numerically unstable.
MORPHOLOGY_FLAGS = POPULATION_FLAGS + ("intrachan_flag",)
#: CHIME sees the sky north of this declination (a configurable choice in the
#: CHIME/FRB system; the catalogue's southernmost burst is at -9.3 deg).
CHIME_DEC_LIMIT_DEG = -11.0
#: CHIME's band edges as they appear in the low_freq/high_freq columns. A burst
#: whose FWTM emission reaches both is "band-filling": its true bandwidth is
#: unknown, only bounded below by the 400 MHz CHIME can see.
BAND_LOW_MHZ = 400.2
BAND_HIGH_MHZ = 800.2
BAND_EDGE_TOL_MHZ = 1.0


# ---------------------------------------------------------------------------
# Dispersion measure bookkeeping
# ---------------------------------------------------------------------------

#: The Milky Way's hot gaseous halo adds DM on top of the disk that NE2001
#: models. Estimates run from ~30 to ~80 pc cm^-3 (Prochaska & Zheng 2019;
#: Cook et al. 2023 from CHIME's lowest-DM bursts). 50 is the value Macquart
#: et al. (2020) adopted; the other two bracket it as a systematic check.
DM_HALO = 50.0
DM_HALO_RANGE = (30.0, 80.0)
#: Upper end of the halo estimate, used when asking whether an observed DM is
#: more than the Milky Way could possibly supply.
DM_HALO_MAX = 80.0
#: Sightline-to-sightline uncertainty in the Milky Way term we subtract, treated
#: as Gaussian in the Macquart fit. NE2001's errors on individual high-latitude
#: sightlines are typically 20-30% of its ~30-60 pc cm^-3, and the halo's column
#: varies with direction by ~10; 20 pc cm^-3 covers both. (The *global* halo
#: level is a separate systematic, bracketed by DM_HALO_RANGE.)
DM_MW_SCATTER = 20.0


# ---------------------------------------------------------------------------
# Macquart relation (DM vs redshift for localized FRBs)
# ---------------------------------------------------------------------------

#: Minimum PATH host-association probability P(O|x) -- the probability that
#: the chosen galaxy, and not a chance neighbour on the sky, is the host. 0.9 is
#: the threshold commonly used to call an association secure.
P_OX_MIN = 0.9
#: Below z ~ 0.01 a host's redshift is dominated by its own peculiar motion
#: rather than cosmic expansion, and the cosmic DM (< 10 pc cm^-3) is smaller
#: than the uncertainty on the Milky Way halo alone. Such bursts carry no
#: information on intergalactic baryons (FRB 20200120E, in a globular cluster of
#: M81 at 3.6 Mpc, is the canonical example).
Z_MIN = 0.01
#: Primordial helium mass fraction (Planck 2018 BBN). With hydrogen and helium
#: both fully ionized below z ~ 3, each baryon mass m_p supplies
#: (1 - Y) + Y/2 = 1 - Y/2 free electrons.
Y_HE = 0.2454
#: Fraction of cosmic baryons in *diffuse ionized* gas (intergalactic medium
#: plus galaxy halos) rather than in stars, stellar remnants and cold ISM.
#: Macquart et al. (2020) adopted 0.844 from the stellar and ISM mass census.
F_DIFFUSE_FIDUCIAL = 0.844
#: Shape indices of the Macquart et al. (2020) DM_cosmic distribution, set by the
#: inner density profile of gas in halos. Both 3 in the original paper.
MACQUART_ALPHA = 3.0
MACQUART_BETA = 3.0

#: Grids for the likelihood fit, flat priors over each range. f_d sweeps well
#: beyond the physically allowed range (it cannot exceed 1) so that the data,
#: not the prior edge, set the answer.
FIT_GRID = {
    "f_d": (0.2, 1.6, 71),
    # log10 of the fluctuation parameter F, where sigma_DM = F / sqrt(z).
    # Macquart et al. (2020) allowed 0.09 <= F <= 0.32, the range of
    # semi-analytic halo-gas models; Baptista et al. (2024) measured
    # log10 F > -0.9 (3 sigma) and found low-z data scatter more than
    # simulations predict. We span 0.04 to 3.2 and let the data choose.
    "log10_F": (-1.4, 0.5, 20),
    # Host DM is lognormal: median exp(mu) from 20 to 800 pc cm^-3 ...
    "mu": (3.0, 6.7, 38),
    # ... with log-width sigma_host.
    "sigma_host": (0.2, 2.0, 19),
}
#: Macquart et al. (2020)'s own, narrower prior on F, for comparison.
MACQUART_2020_LOG10_F = (-1.0458, -0.4949, 8)  # F from 0.09 to 0.32
#: Points in the numerical convolution over the host's DM contribution.
FIT_N_CONV = 400
#: Gauss-Hermite nodes for integrating over the Milky Way scatter.
FIT_N_MW_NODES = 7

#: FRB 20190520B lives behind two foreground galaxy clusters and has a host DM
#: of ~900 pc cm^-3 (Niu et al. 2022; Lee et al. 2023), far above the rest. The
#: model's lognormal tail can absorb it, but we also refit without it.
KNOWN_HOST_OUTLIERS = ("FRB20190520B",)
