"""Project configuration: the observing epoch, data sources, and every tunable choice.

"Today" in this project is a pinned moment, ``DEFAULT_EPOCH``, so that the numbers
quoted in the README can be regenerated exactly. Every run can also be pointed at
another moment (``python run.py --epoch now`` for the Sun as it is right now); each
epoch gets its own cache under ``data/<tag>/`` and, unless it is the default, its
own figure folder ``figures/<tag>/``.

Everything the analysis depends on is set here, with the reason for each choice.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_ROOT = PROJECT_DIR / "data"
FIGURE_ROOT = PROJECT_DIR / "figures"

# --- The epoch: what "today" means ---------------------------------------------------

#: 2026-09-29 20:00 UTC, the latest hour with a full set of SDO images when this
#: project was built. SDO's near-real-time pipeline delivers images 15-20 minutes
#: after they are taken, so ``--epoch now`` resolves to the last complete hour
#: before that.
DEFAULT_EPOCH = "2026-09-29T20:00"

#: The week of history behind the snapshot: X-ray flux, flares, and the images
#: from which we measure the Sun's rotation.
WEEK_DAYS = 7


def parse_epoch(text: str | None) -> datetime:
    """Turn ``None``/``"now"``/an ISO time into a UTC datetime on the hour."""
    if text is None:
        text = DEFAULT_EPOCH
    if text == "now":
        t = datetime.now(timezone.utc) - timedelta(minutes=40)
    else:
        t = datetime.fromisoformat(text)
        t = t.replace(tzinfo=timezone.utc) if t.tzinfo is None else t.astimezone(timezone.utc)
    return t.replace(minute=0, second=0, microsecond=0)


@dataclass(frozen=True)
class EpochPaths:
    """Where one epoch's data and figures live."""

    epoch: datetime
    tag: str
    data: Path
    figures: Path

    @property
    def snapshot(self) -> Path:
        return self.data / "snapshot"

    @property
    def week_aia(self) -> Path:
        return self.data / "week" / f"aia{WEEK_AIA_CHANNEL}"

    @property
    def week_hmi(self) -> Path:
        return self.data / "week" / "hmi"

    @property
    def continuum(self) -> Path:
        return self.data / "hmi_continuum"

    @property
    def goes(self) -> Path:
        return self.data / "goes"

    @property
    def flare(self) -> Path:
        return self.data / "flare"

    @property
    def manifest(self) -> Path:
        return self.data / "manifest.json"

    @property
    def results(self) -> Path:
        return self.data / "results.json"

    @property
    def cache(self) -> Path:
        """Intermediate products (reprojected maps, tracked regions), all rebuildable."""
        return self.data / "derived"

    @property
    def week_start(self) -> datetime:
        return self.epoch - timedelta(days=WEEK_DAYS)

    def ensure(self) -> "EpochPaths":
        for d in (self.data, self.figures):
            d.mkdir(parents=True, exist_ok=True)
        return self


def paths_for(epoch: datetime) -> EpochPaths:
    tag = epoch.strftime("%Y%m%dT%H%M")
    is_default = epoch == parse_epoch(DEFAULT_EPOCH)
    return EpochPaths(
        epoch=epoch,
        tag=tag,
        data=DATA_ROOT / tag,
        figures=FIGURE_ROOT if is_default else FIGURE_ROOT / tag,
    )


# --- Data sources -----------------------------------------------------------------------
#
# SDO's near-real-time "synoptic" images: every AIA channel every 3 minutes, taken
# from the 4096x4096 level-1 frames, registered to level 1.5 (Sun centred, solar
# north up, 0.6"/pixel) and averaged 4x4 down to 1024x1024 pixels of 2.4". Pixel
# values are DN per original 0.6" pixel, NOT per second: we divide by EXPTIME.
# Latency 15-20 minutes; files persist, so any past moment since 2023 works too.
AIA_NRT_ROOT = "https://jsoc1.stanford.edu/data/aia/synoptic/nrt"

#: HMI line-of-sight magnetograms, 720-s averages, also binned to 1024x1024
#: (2.0"/pixel), one per hour. Recent hours are near-real-time ("_nrt") files; after
#: a few days they are replaced by the definitive version of the same hour.
HMI_FITS_ROOT = "http://jsoc.stanford.edu/data/hmi/fits"

#: NOAA SWPC's operational GOES X-ray feeds: 1-minute fluxes from the primary
#: satellite and the event list SWPC derives from them. These are rolling 7-day
#: windows, so they only cover an epoch within the last week of the download.
SWPC_XRAYS_URL = "https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json"
SWPC_FLARES_URL = "https://services.swpc.noaa.gov/json/goes/primary/xray-flares-7-day.json"
SWPC_REGIONS_URL = "https://services.swpc.noaa.gov/json/solar_regions.json"

#: NOAA NCEI's science-quality GOES-R XRS archive (1-minute averages, one netCDF
#: file per day, ~1 day latency). Used for epochs the SWPC window no longer covers.
NCEI_XRS_ROOT = (
    "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/"
    "goes{sat}/l2/data/xrsf-l2-avg1m_science"
)
#: NCEI's science flare list ("flsum"): start, peak, end and class of each flare.
NCEI_FLSUM_ROOT = NCEI_XRS_ROOT.replace("avg1m", "flsum")

#: CHIANTI-based GOES/XRS temperature response table (White, Thomas & Schwartz
#: 2005, updated for GOES-R), as distributed with SolarSoft and read by
#: sunkit-instruments. Maps the ratio of the two X-ray channels to a temperature.
GOES_RESPONSE_URLS = (
    "https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits",
    "https://sohoftp.nascom.nasa.gov/solarsoft/gen/idl/synoptic/goes/goes_chianti_response_latest.fits",
)

# --- The snapshot: the Sun at the epoch, in every channel ------------------------------

#: How far from the epoch a snapshot frame may be. AIA's 4500 channel is only
#: taken on the hour; everything else every 3 minutes.
SNAPSHOT_TOLERANCE_MIN = 6
HMI_TOLERANCE_MIN = 90


@dataclass(frozen=True)
class Channel:
    """One of AIA's wavelengths and what it sees (Lemen et al. 2012, Table 1)."""

    wavelength: int  # angstrom
    ions: str
    log_t: tuple  # characteristic temperature(s), log10 K
    region: str  # where in the atmosphere the emission forms
    cmap: str  # sunpy's standard SDO colour table
    stretch: str  # display stretch: "log" or "asinh"


#: All ten AIA channels, ordered by the temperature of the gas that dominates each
#: one on the quiet Sun: a walk upward from the visible surface to flaring loops.
#: (131 and 193 each contain a second, much hotter ion that lights up in flares.)
CHANNELS = (
    Channel(4500, "continuum", (3.7,), "photosphere", "sdoaia4500", "asinh"),
    Channel(1700, "continuum", (3.7,), "temperature minimum, photosphere", "sdoaia1700", "asinh"),
    Channel(1600, "C IV + continuum", (5.0, 3.7), "transition region, upper photosphere",
            "sdoaia1600", "asinh"),
    Channel(304, "He II", (4.7,), "chromosphere, transition region", "sdoaia304", "log"),
    Channel(131, "Fe VIII, Fe XXI", (5.6, 7.0), "transition region, flaring corona",
            "sdoaia131", "log"),
    Channel(171, "Fe IX", (5.8,), "quiet corona, upper transition region", "sdoaia171", "log"),
    Channel(193, "Fe XII, Fe XXIV", (6.2, 7.3), "corona, hot flare plasma", "sdoaia193", "log"),
    Channel(211, "Fe XIV", (6.3,), "active-region corona", "sdoaia211", "log"),
    Channel(335, "Fe XVI", (6.4,), "active-region corona", "sdoaia335", "log"),
    Channel(94, "Fe XVIII", (6.8,), "flaring corona", "sdoaia94", "log"),
)
CHANNEL_BY_WAVELENGTH = {c.wavelength: c for c in CHANNELS}

#: The AIA quality word: bit 30 only marks near-real-time ("quicklook") data, which
#: is what we use. Any other bit set means a known problem with the frame.
AIA_QUALITY_ALLOWED_MASK = 1 << 30

# --- The week: rotation and the time-lapse ---------------------------------------------

#: One magnetogram and one AIA image every 3 hours for the week: 56 of each. The
#: Sun turns 1.6 degrees in 3 hours, which keeps the time-lapse smooth, and 24-hour
#: pairs for the rotation measurement overlap at 8 different phases of the day.
WEEK_CADENCE_HOURS = 3
#: 193 A shows active regions, coronal holes and loops in a single channel.
WEEK_AIA_CHANNEL = 193

# --- Differential rotation from magnetogram cross-correlation ---------------------------
#
# The classic magnetic-feature method (Snodgrass 1983; Komm, Howard & Harvey 1993):
# remap each magnetogram onto a Carrington longitude-latitude grid, then find how far
# the pattern in each latitude band has drifted after a day.

#: Carrington grid resolution. An HMI 1k pixel is 2.0" = 1.45 Mm = 0.12 degrees of
#: heliographic arc at disc centre, so 0.2 degrees keeps the full detail there.
ROT_GRID_DEG = 0.2
#: Only use the disc within this many degrees of its centre: nearer the limb the
#: line-of-sight field of vertical flux tubes fades as cos(angle), and pixels grow.
ROT_MAX_DISK_ANGLE_DEG = 60.0
#: Time between the two magnetograms of a pair. The magnetic network is recycled
#: within a day or two, so the pattern decorrelates: this week, pairs 12 h apart
#: kept a correlation of 0.3-0.4 at every latitude to +-47 deg, while 24-h pairs
#: fell to ~0.25 and lost everything beyond +-25 deg; 48-h pairs kept nothing. The
#: drift to measure is small either way (the Carrington grid co-rotates), and the
#: other separations are run too, as a check on the method (``ROT_ALT_SEPARATIONS_H``).
ROT_PAIR_SEPARATION_H = 12
ROT_ALT_SEPARATIONS_H = (6, 18, 24)
#: Latitude bands for the rotation profile.
ROT_BAND_EDGES_DEG = tuple(range(-50, 51, 5))
#: Active regions are masked, as Komm et al. (1993) did, because a large region
#: would dominate its band's correlation and regions rotate differently from the
#: small-scale network. Mask = smoothed |B_r| above this many gauss, grown by a
#: margin so the region's surrounding plage goes too.
ROT_AR_MASK_GAUSS = 60.0
ROT_AR_MASK_SMOOTH_PX = 4.0
ROT_AR_MASK_GROW_DEG = 2.0
#: Clip the network field so a few strong elements cannot dominate a band.
ROT_FIELD_CLIP_GAUSS = 200.0
#: Largest drift searched, in degrees of Carrington longitude per pair. The
#: Carrington frame itself turns at 14.18 deg/day, so real drifts are tiny.
ROT_MAX_SHIFT_DEG = 4.0

#: The Carrington frame's sidereal rotation rate: 360 degrees in 25.38 days, as
#: adopted by the IAU (W = 84.176 + 14.1844 d).
CARRINGTON_RATE_DEG_PER_DAY = 14.1844

# --- Active regions on the week's magnetograms ------------------------------------------

#: A region is a patch where the smoothed |B_los| exceeds this, larger than
#: ``AR_MIN_AREA_MSH`` millionths of a hemisphere (the unit NOAA uses for spot areas).
AR_DETECT_GAUSS = 50.0
AR_DETECT_SMOOTH_PX = 4.0
AR_MIN_AREA_MSH = 150.0
#: A region's two polarities are often separated by weak field; patches closer than
#: this (pixels of 2", so 24" = 17 Mm) are grouped into one region before measuring.
AR_GROUP_DILATE_PX = 12
#: Detections on successive magnetograms belong to the same region if, after
#: allowing for differential rotation, their centres lie within this distance.
AR_LINK_DEG = 5.0
#: A region needs at least this many hours of tracking for a rotation rate.
AR_MIN_TRACK_H = 48.0
#: NOAA region numbers are attached when a SWPC position lies within this distance.
AR_NOAA_MATCH_DEG = 8.0

# --- Limb darkening (HMI continuum) ------------------------------------------------------

#: HMI's continuum images are not in the near-real-time FITS feed, so we take the
#: most recent definitive one via the Virtual Solar Observatory: the photosphere's
#: limb darkening does not change from one day to the next.
CONTINUUM_SEARCH_BACK_DAYS = 12
#: Fit range in mu = cos(angle from disc centre). Below mu = 0.1 (the outer 0.5% of
#: the radius) the image's blur and scattered light dominate the true profile.
LIMB_MU_MIN = 0.1
#: Radial bins for the profile, in mu.
LIMB_N_BINS = 60
#: HMI observes the Fe I 6173.3 A line; the continuum it reports lies alongside.
HMI_WAVELENGTH_A = 6173.3

# --- The corona above the limb -----------------------------------------------------------

#: Heights above the photospheric limb used to measure the coronal scale height.
#: Below ~10 Mm spicules and the transition region dominate; above ~80 Mm the
#: signal fades toward instrumental scattered light.
CORONA_FIT_MIN_MM = 12.0
CORONA_FIT_MAX_MM = 80.0
#: Position-angle sectors around the limb, each fitted separately.
CORONA_N_SECTORS = 36
#: The channels whose off-limb emission we fit: the four iron-line coronal
#: channels with enough signal above the quiet limb.
CORONA_CHANNELS = (171, 193, 211, 335)

# --- Flares --------------------------------------------------------------------------------

#: GOES flare classes by peak 1-8 A flux (W/m^2).
FLARE_CLASSES = (("A", 1e-8), ("B", 1e-7), ("C", 1e-6), ("M", 1e-5), ("X", 1e-4))
#: The flare studied in detail: the one of at least this class with the largest
#: soft X-ray energy (background-subtracted 1-8 A fluence) in the week.
FLARE_MIN_CLASS_FLUX = 1e-6
#: Background for a flare: the median 1-8 A flux in this window before it starts.
FLARE_BACKGROUND_MIN = 10
#: AIA frames around the flare: from this long before the GOES start to this long
#: after its end. Flare loops keep cooling through the EUV channels for hours after
#: the X-rays have faded: the 2 MK channels were still brightening 60 min after the
#: end of this week's flare, so the window runs to 3 hours.
FLARE_BEFORE_MIN = 15
FLARE_AFTER_MIN = 180
#: The channels followed through the flare, from the chromosphere to 10 MK.
FLARE_CHANNELS = (1600, 304, 131, 94, 335, 211, 193, 171)
#: Cut-out size around the flare site. Full frames are downloaded one at a time,
#: cut, and discarded, so the cache holds ~10 MB instead of ~270 MB.
FLARE_CUTOUT_HALF_WIDTH_ARCSEC = 150.0
#: Light curves are summed over this radius around the flare site.
FLARE_APERTURE_ARCSEC = 45.0

# --- Physical constants and published comparisons ------------------------------------

#: IAU 2015 nominal solar values (Prsa et al. 2016).
R_SUN_M = 6.957e8
T_EFF_K = 5772.0
GM_SUN = 1.3271244e20  # m^3 s^-2
#: SDO's own FITS headers use this radius (RSUN_REF) for the solar limb.
R_SUN_REF_M = 6.96e8

#: Helium-to-hydrogen number ratio of the corona: sets the mean mass per particle.
CORONA_HE_ABUNDANCE = 0.085

#: Published differential rotation laws, omega = A + B sin^2(lat) + C sin^4(lat) in
#: sidereal degrees/day, as compiled by Sudar et al. (2015, A&A 575, A63, Table 2).
PUBLISHED_ROTATION = {
    "snodgrass1983": dict(label="Snodgrass (1983), magnetic features", tracer="magnetic",
                          A=14.307, B=-1.98, C=-2.15, sA=0.005, sB=0.06, sC=0.11),
    "komm1993": dict(label="Komm et al. (1993), magnetic features", tracer="magnetic",
                     A=14.42, B=-2.00, C=-2.09, sA=0.02, sB=0.13, sC=0.15),
    "snodgrass1984": dict(label="Snodgrass (1984), Doppler (surface gas)", tracer="doppler",
                          A=14.05, B=-1.49, C=-2.61),
    "balthasar1986": dict(label="Balthasar et al. (1986), sunspot groups", tracer="sunspots",
                          A=14.551, B=-2.87, C=0.0, sA=0.006, sB=0.06),
    "schou1998": dict(label="Schou et al. (1998), helioseismology", tracer="interior",
                      A=14.16, B=-1.63, C=-2.52),
}

#: Hale's polarity law for Cycle 25: leading (westward) spots are positive in the
#: north and negative in the south, the reverse of Cycle 24.
CYCLE25_LEADING_POLARITY = {"N": +1, "S": -1}
