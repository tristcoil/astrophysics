"""Project configuration: paths, data sources, conventions and every sample cut.

L2 is a census of exploding stars. Its backbone is the Zwicky Transient
Facility (ZTF) Bright Transient Survey (BTS): since 2018 the BTS team has taken
a spectrum of nearly every transient in the northern sky that peaked brighter
than 18.5 mag, so its catalogue is a magnitude-limited, spectroscopically
complete sample of ~10,000 supernovae. From it we measure how common each kind
of supernova is, how bright, and in which galaxies. Light curves of the Type II
supernovae (ALeRCE) give their plateaus, and published pre-explosion images of
nearby supernovae give the stars that blew up, and the stars that apparently
did not.

Every cut lives here, next to the reason for it, so any number in the README
can be traced back to exactly which objects it came from.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
#: The BTS Sample Explorer's CSV exports, kept exactly as served.
BTS_DIR = DATA_DIR / "bts"
#: One JSON file per supernova: ZTF alert photometry from the ALeRCE broker.
ALERCE_DIR = DATA_DIR / "alerce"
#: SDSS galaxies around each supernova (one CSV per batch of positions).
SDSS_HOST_DIR = DATA_DIR / "sdss_hosts"
#: MPA-JHU masses and star-formation rates of all SDSS galaxies at z < 0.1.
SDSS_BACKGROUND_DIR = DATA_DIR / "sdss_background"
#: Cosmicflows-4 (Tully et al. 2023): redshift-independent distances.
CF4_DIR = DATA_DIR / "cf4"
#: Galactic dust map (Schlegel, Finkbeiner & Davis 1998), HEALPix version.
DUST_DIR = DATA_DIR / "dust"
#: LaTeX sources of the papers whose tables we parse ourselves.
PAPERS_DIR = DATA_DIR / "papers"
#: ZTF public-survey exposure history from IRSA (optional, see rates.py).
ZTF_EXPOSURE_DIR = DATA_DIR / "ztf_exposures"
#: Intermediate products of the analysis (host matches, light-curve fits).
DERIVED_DIR = DATA_DIR / "derived"
RESULTS_PATH = DATA_DIR / "results.json"
FIGURE_DIR = PROJECT_DIR / "figures"

for _d in (DATA_DIR, BTS_DIR, ALERCE_DIR, SDSS_HOST_DIR, SDSS_BACKGROUND_DIR, CF4_DIR,
           DUST_DIR, PAPERS_DIR, ZTF_EXPOSURE_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Where the data comes from
# ---------------------------------------------------------------------------

#: The BTS Sample Explorer (Perley et al. 2020). It is regenerated from the
#: survey database every few hours; we keep the snapshot we analysed.
BTS_EXPLORER_URL = "https://sites.astro.caltech.edu/ztf/bts/explorer.php"
#: Query for *every* saved candidate (no quality, purity or magnitude cut, all
#: classes, host photometry and sample flags shown), so every cut is ours to
#: apply and document. Field names are the explorer's own form fields.
BTS_QUERY_ALL = {
    "f": "s", "subsample": "all", "classstring": "", "classexclude": "",
    "hostphot": "y", "samp": "y", "discname": "y", "showdiscdate": "y",
    "ztflink": "alerce", "lastdet": "",
    **{f"{side}{key}": "" for side in ("start", "end")
       for key in ("savedate", "peakdate", "ra", "dec", "z", "dur", "rise", "fade", "peakmag",
                   "absmag", "hostabs", "hostcol", "b", "av")},
    "sort": "id",
    **{f"show{key}": "y" for key in ("savedate", "peakdate", "lastdate", "ra", "dec", "z", "rise",
                                     "fade", "peakmag", "lastmag", "absmag")},
    "format": "csv",
}
#: The same query restricted to SDSS photometric fields ("Require SDSS field").
#: Perley et al. trust the explorer's host photometry only there: Pan-STARRS
#: catalogue photometry shreds nearby galaxies and is not aperture-matched.
BTS_QUERY_SDSS = {**BTS_QUERY_ALL, "sdss": "y"}

#: ALeRCE broker (Förster et al. 2021): every public ZTF alert, re-served as a
#: light curve per object. No account needed.
ALERCE_LC_URL = "https://api.alerce.online/ztf/v1/objects/{oid}/lightcurve"

#: SkyServer's public SQL endpoint for SDSS Data Release 17, which still
#: carries the MPA-JHU value-added tables (galSpec*).
SKYSERVER_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
SKYSERVER_CUTOUT_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/ImgCutout/getjpeg"

#: Cosmicflows-4 on the CDS FTP server (VizieR J/ApJ/944/94). A browser-like
#: user agent gets an anti-bot challenge there; a plain one gets the files.
CF4_BASE_URL = "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/944/94/"
CF4_FILES = ("ReadMe", "table2.dat.gz", "table4.dat.gz")

#: SFD (1998) E(B-V) as a HEALPix map (NASA LAMBDA), used for the survey's
#: average loss of volume to Galactic dust.
SFD_URL = "https://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits"

#: IRSA's table of ZTF science exposures. Optional: when the service answers we
#: measure the survey's active sky area ourselves (see rates.py); when it does
#: not, the rates fall back on Perley et al.'s published coverage.
IRSA_TAP_URL = "https://irsa.ipac.caltech.edu/TAP/sync"
#: ZTF field centres (the survey's fixed grid), from the ZTF team's GitHub.
ZTF_FIELDS_URL = ("https://raw.githubusercontent.com/ZwickyTransientFacility/"
                  "ztf_information/master/field_grid/ZTF_Fields.txt")

#: arXiv LaTeX sources whose tables we parse (tarball URL = arxiv.org/src/<id>).
ARXIV_SOURCES = {
    # Smartt 2015, PASA 32, e016: the compilation of progenitor detections and
    # limits behind "13 expected, about one found".
    "smartt2015": "1504.02635",
    # Davies & Beasor 2018, MNRAS 474, 2116: the same supernovae with
    # bolometric corrections measured on red supergiants near death.
    "davies_beasor2018": "1709.06116",
}


# ---------------------------------------------------------------------------
# The BTS catalogue: flags and classes
# ---------------------------------------------------------------------------

#: The explorer's ``qualcode`` string, one character per quality condition
#: (decoded by querying the explorer with each condition alone and matching
#: the returned objects to the bit, 100% agreement for every bit):
QUALCODE_BITS = ("coverage", "visibility", "alert_filter", "clean_reference",
                 "peak_after_2018may", "low_extinction")
#: ``purecode``: SN-like light curve, galaxy cross-match (either one passes).
PURECODE_BITS = ("sn_like_lightcurve", "galaxy_match")

#: How the TNS classifications are grouped. Following Perley et al. (2020,
#: their Fig. 7): thermonuclear; hydrogen-rich (ordinary II and IIb);
#: interacting (IIn, SLSN-II, Ibn/Icn, Ia-CSM); hydrogen-poor (Ib/c, Ic-BL,
#: SLSN-I). We keep a finer split for the luminosity functions.
CLASS_MAP = {
    "SN Ia": "Ia", "SN Ia-91T": "Ia", "SN Ia-91bg": "Ia", "SN Ia-pec": "Ia",
    "SN Iax": "Ia", "SN Ia-SC": "Ia",
    "SN II": "II", "SN IIP": "II", "SN IIL": "II", "SN II-pec": "II",
    "SN IIb": "IIb",
    "SN IIn": "IIn", "SN IIn-pec": "IIn",
    "SN Ib": "Ibc", "SN Ic": "Ibc", "SN Ib/c": "Ibc", "SN Ib-pec": "Ibc", "SN Ic-pec": "Ibc",
    "SN Ic-BL": "IcBL",
    "SN Ibn": "Ibn/Icn", "SN Icn": "Ibn/Icn",
    "SN Ia-CSM": "Ia-CSM",
    "SLSN-I": "SLSN-I", "SLSN-II": "SLSN-II",
    "TDE": "TDE", "TDE-He": "TDE", "TDE-H-He": "TDE", "TDE-featureless": "TDE",
}
#: Fine class -> broad family used for fractions and rates.
FAMILY = {
    "Ia": "Ia", "II": "II", "IIb": "II", "IIn": "IIn", "Ibc": "Ibc", "IcBL": "Ibc",
    "Ibn/Icn": "IIn", "Ia-CSM": "IIn", "SLSN-I": "SLSN", "SLSN-II": "SLSN", "TDE": "TDE",
}
#: Which families count as core collapse. Perley et al. count superluminous
#: supernovae as core collapse; Ia-CSM (a thermonuclear explosion inside dense
#: gas) is not, but there are so few that the choice is invisible.
CORE_COLLAPSE_CLASSES = ("II", "IIb", "IIn", "Ibc", "IcBL", "Ibn/Icn", "SLSN-I", "SLSN-II")
#: The display order of the families (fixed, so a colour always means one family).
FAMILY_ORDER = ("Ia", "II", "Ibc", "IIn", "SLSN", "TDE")
FAMILY_LABEL = {
    "Ia": "Type Ia (thermonuclear)", "II": "Type II (+ IIb)", "Ibc": "Type Ib/c (+ Ic-BL)",
    "IIn": "Interacting (IIn, Ibn, Icn, Ia-CSM)", "SLSN": "Superluminous", "TDE": "Tidal disruption",
}


# ---------------------------------------------------------------------------
# Conventions inherited from the BTS (so our numbers line up with theirs)
# ---------------------------------------------------------------------------

#: Perley et al. (2020) compute absolute magnitudes with H0 = 70, Om = 0.3,
#: OL = 0.7, from heliocentric redshifts, with a "K-correction" of
#: +2.5 log10(1+z) (the bandwidth term for a flat spectrum). We reproduce the
#: explorer's own M_abs column to 0.003 mag with these choices (see bts.py).
H0 = 70.0
OMEGA_M = 0.3
#: Galactic extinction in the ZTF bands per unit A_V, fitted to reproduce the
#: explorer's absolute magnitudes (median ratio over 4,300 events with
#: A_V > 0.15). Schlafly & Finkbeiner (2011) give 1.205 and 0.833 for SDSS g, r.
A_G_PER_AV = 1.208
A_R_PER_AV = 0.837

#: Below this heliocentric velocity a redshift is a poor distance: peculiar
#: motions of ~300 km/s are a large fraction of cz. SN 2023ixf in M101 is the
#: extreme case: z = 0.0008 puts it at 3.4 Mpc instead of the Cepheid 6.85 Mpc,
#: 1.5 mag too faint. Hosts below this velocity get a Cosmicflows-4 distance.
CF4_MAX_CZ_KMS = 4000.0
#: Search radius for the Cosmicflows-4 host: projected kpc at the SN's own
#: redshift distance (floored at 2'), and the allowed velocity mismatch.
CF4_MATCH_KPC = 40.0
CF4_MATCH_MIN_ARCMIN = 2.0
CF4_MATCH_DV_KMS = 600.0
#: Cosmicflows-4 distances sit on the Cepheid/TRGB scale, which implies
#: H0 = 74.6 (Tully et al. 2023). Everything else here uses H0 = 70, so CF4
#: moduli are shifted by 5 log10(74.6/70) = +0.138 mag to the same scale
#: (rates are quoted per h70^-3 Gpc^3, as the BTS papers quote them).
CF4_H0 = 74.6
#: The CMB dipole (Planck 2018): used to move redshifts into the CMB frame.
CMB_DIPOLE_L_DEG = 264.021
CMB_DIPOLE_B_DEG = 48.253
CMB_DIPOLE_KMS = 369.82
C_KMS = 299792.458


# ---------------------------------------------------------------------------
# The statistical sample
# ---------------------------------------------------------------------------

#: Classification limit of the survey. Perley et al.: 93% of transients that
#: peak brighter than this are classified; we measure the completeness per
#: magnitude bin ourselves.
M_LIM = 18.5
#: Peak dates of the statistical sample. The start is the public start of the
#: BTS (Perley et al.). The end is the last month whose candidates are >= 90%
#: classified: from November 2025 completeness falls (87%, 84%, then 30% in
#: January 2026) as the classification effort wound down.
STAT_PEAK_START = "2018-06-01"
STAT_PEAK_END = "2025-10-31"
#: Perley et al.'s own window, for the replication of their rates.
PERLEY_PEAK_START = "2018-06-01"
PERLEY_PEAK_END = "2020-07-15"

#: Perley et al. (2020) rate ingredients, used as published for the replication.
PERLEY_F_SKY = 0.35           # 14,400 deg2 average active area
PERLEY_F_EXT = 0.82           # average loss of volume to Galactic dust
PERLEY_F_REC = 0.60           # 412 / 690 easy SNe Ia pass the quality cut
PERLEY_T_YR = 2.12
PERLEY_FCL_BRIGHT_MAG = 17.2  # classification efficiency 1.0 brighter than this,
PERLEY_FCL_AT_LIMIT = 0.9     # falling linearly to 0.9 at m = 18.5
#: Perley et al.'s "easy" SNe Ia, whose recovery measures f_rec: peak absolute
#: magnitude brighter than -18.5 and apparent peak brighter than 18.0.
EASY_IA_MABS = -18.5
EASY_IA_MAPP = 18.0
#: Integration limits of the volumetric rates: SNe Ia brighter than -16.5,
#: core-collapse brighter than -14, superluminous brighter than -21
#: (Perley et al.'s choices).
RATE_MLIM_IA = -16.5
RATE_MLIM_CC = -14.0
RATE_MLIM_SLSN = -21.0

#: The BTS footprint as applied by its own cuts: Galactic latitude, the
#: declination reached by the visibility cut, and the extinction cut.
FOOTPRINT_MIN_ABS_B = 7.0
FOOTPRINT_MIN_DEC = -26.0
FOOTPRINT_MAX_AV = 1.0
#: Pan-STARRS covers north of -30 deg; the visibility cut (30 deg elevation for
#: two hours, 30 days after peak) removes everything south of -26 deg.

#: Bootstrap resamples for every rate and fraction.
N_BOOT = 2000
RNG_SEED = 20261004


# ---------------------------------------------------------------------------
# Host galaxies
# ---------------------------------------------------------------------------

#: Search radius around each supernova: the smaller of 90" and 30 projected kpc
#: at its redshift (the BTS team's own limits), floored at 20".
HOST_SEARCH_MAX_ARCSEC = 90.0
HOST_SEARCH_MAX_KPC = 30.0
HOST_SEARCH_MIN_ARCSEC = 20.0
#: Host = the galaxy with the smallest "directional light radius" distance
#: (separation in units of the galaxy's elliptical radius toward the SN;
#: Sullivan et al. 2006, Gupta et al. 2016). Beyond 4 DLR there is no host.
HOST_MAX_DLR = 4.0
#: Faint limit for candidate host galaxies in SDSS photometry.
HOST_MAX_RMAG = 22.5
#: A host with a spectrum must agree with the SN redshift: SN-template
#: redshifts carry ~0.005 errors (Fremling et al. 2020).
HOST_DZ_MAX = 0.01
#: Batch size for the SkyServer cross-match (URL length limit ~100 positions).
SKYSERVER_BATCH = 80
#: Redshift range for the SFR-mass background and the host analysis: the SDSS
#: Main Galaxy Sample is reasonably complete for log M > 9.5 here.
HOST_Z_MIN = 0.015
HOST_Z_MAX = 0.1
#: MPA-JHU spectra live on SDSS-I/II plates 266-2974; fetched 100 plates at a time.
GALSPEC_PLATE_MIN = 266
GALSPEC_PLATE_MAX = 2974
GALSPEC_PLATE_CHUNK = 100
#: Quenched = specific SFR below this (log yr^-1), the usual dividing line
#: ~1 dex below the star-forming main sequence at z ~ 0.05.
LOG_SSFR_QUENCHED = -11.0


# ---------------------------------------------------------------------------
# Type II light curves
# ---------------------------------------------------------------------------

#: Which objects get a light curve from ALeRCE: every Type II (II, IIP, IIL,
#: II-pec) for the plateau analysis, plus the brightest examples of every
#: class for the gallery and the explorer.
LC_TYPE_II_CLASSES = ("SN II", "SN IIP", "SN IIL", "SN II-pec")
LC_GALLERY_PER_CLASS = 16
#: Polite parallelism for the broker.
ALERCE_WORKERS = 3

#: Explosion epoch: halfway between the last non-detection and the first
#: detection, if they are at most this many rest-frame days apart and the
#: limit was deeper than the first detection.
EXPLOSION_MAX_GAP_DAYS = 10.0
#: Lyman, Bersier & James (2014) bolometric correction for Type II supernovae:
#: BC_g = c0 + c1 (g-r) + c2 (g-r)^2, valid for -0.2 < g-r < 1.3 (rms 0.036 mag).
LYMAN_BC_II = (0.053, -0.089, -0.736)
LYMAN_BC_COLOR_RANGE = (-0.2, 1.3)
M_BOL_SUN = 4.74
#: Progenitor radius assumed when turning plateaus into ejecta masses
#: (Goldberg et al. 2019 inversions): red supergiants of log L ~ 4.5-5.2 and
#: T ~ 3,500-4,000 K have R ~ 400-1,000 R_sun.
#: Plateaus brighter than this at day 50 are left out of the energy/mass
#: inversion: no recombination-powered plateau is that luminous (Goldberg et
#: al.'s models reach ~10^42.6 erg/s), so extra power from collision with
#: surrounding gas is the likely source (one event, SN 2019lkw at -19.2).
LC_MAX_PLATEAU_MR = -18.5
R_PROGENITOR_RSUN = 500.0
R_PROGENITOR_RANGE_RSUN = (300.0, 1000.0)

# ---------------------------------------------------------------------------
# Progenitors
# ---------------------------------------------------------------------------

#: Salpeter (1955) initial mass function slope, dN/dM ~ M^-2.35.
SALPETER_SLOPE = -2.35
#: Smartt (2015)'s mass range for the "13 expected" argument.
SMARTT_M_MIN = 8.0
SMARTT_M_MAX = 100.0
SMARTT_M_SPLIT = 18.0
#: Empirical upper luminosity of red supergiants (Humphreys-Davidson limit for
#: cool supergiants; Davies, Crowther & Beasor 2018).
LOGL_HD_LIMIT = 5.5
#: Davies & Beasor (2018) mass-luminosity relation of the STARS models at the
#: start of neon burning: log M = a + b log L (valid 8-25 M_sun).
DB18_MLR = (-1.2417, 0.47868)
#: Beasor, Smith & Jencson (2024), testing the progenitor recipe on 139 red
#: supergiants in M31 whose full-SED luminosities are known: Smartt's assumed
#: early-M bolometric correction (BC_I = 0.9) underestimates the luminosity by
#: 0.3 dex on average. A late-M correction (BC_I = 0, as Davies & Beasor 2018
#: adopted) removes the mean offset (-0.08 dex) but leaves a 0.2 dex scatter,
#: which they add to every single-band (I/F814W) progenitor luminosity.
BEASOR24_OFFSET_DEX = 0.3
BEASOR24_EXTRA_ERR_DEX = 0.2
N_MC_PROGENITORS = 20000
