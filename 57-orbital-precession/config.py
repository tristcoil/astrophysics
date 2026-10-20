"""Project configuration: paths, physical conventions, the N-body experiments, and
the published numbers we grade ourselves against.

The backbone is JPL's planetary ephemeris DE440 (Park et al. 2021, AJ 161, 105):
the same numerical integration of the Solar System that JPL navigates spacecraft
with, fitted to a century of radar, spacecraft and lunar laser ranging. It supplies

* the **initial conditions** of our own N-body runs (every planet's position and
  velocity at J2000.0),
* the **masses** (GM values) and the Sun's J2 and spin, read from its header file,
* the **answer key**: Mercury's orbit from 1550 to 2650 with all the physics JPL
  includes (general relativity, the Sun's oblateness and spin, 343 asteroids, ...).

Everything else in the project -- which forces to switch on, how to measure a
precession rate -- is ours, and is checked against DE440 and the literature.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
DE440_DIR = RAW_DIR / "de440"
DE440_BSP = DE440_DIR / "de440.bsp"
DE440_HEADER = DE440_DIR / "header.440"
DERIVED_DIR = DATA_DIR / "derived"
RUNS_DIR = DERIVED_DIR / "runs"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, DE440_DIR, DERIVED_DIR, RUNS_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Where the data come from
# --------------------------------------------------------------------------

DE440_BSP_URL = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440.bsp"
DE440_HEADER_URL = "https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440"
HORIZONS_API = "https://ssd.jpl.nasa.gov/api/horizons.api"

# --------------------------------------------------------------------------
# Units and conventions
# --------------------------------------------------------------------------

#: Julian date of the standard epoch J2000.0 (2000 January 1, 12:00 TDB).
J2000_JD = 2451545.0
#: A Julian century, the unit every precession rate in this project is quoted in.
CENTURY_DAYS = 36525.0
YEAR_DAYS = 365.25
ARCSEC_PER_RAD = 180.0 / 3.141592653589793 * 3600.0

#: The "ecliptic of J2000" that JPL Horizons and the astronomical literature use:
#: the ICRF rotated about its x-axis by the IAU 1976 obliquity. Horizons' own
#: header states this value (84381.448 arcsec). It is a *fixed* plane: unlike the
#: ecliptic of date, it does not move as the Earth's orbit is perturbed.
OBLIQUITY_J2000_ARCSEC = 84381.448

#: The Sun's rotation pole in the ICRF (IAU Working Group on Cartographic
#: Coordinates and Rotational Elements; Archinal et al. 2018, Celest. Mech. Dyn.
#: Astron. 130, 22). The axis of the Sun's oblateness (J2) and of its spin.
SUN_POLE_RA_DEG = 286.13
SUN_POLE_DEC_DEG = 63.87

#: The span DE440 covers is 1549-12-21 to 2650-01-25 (JD 2287184.5-2688976.5).
#: Every run starts at J2000.0 and is integrated backwards and forwards to the
#: edges of that span, so the whole of DE440 can be compared epoch by epoch.
SPAN_BACK_YR = 449.5
SPAN_FORWARD_YR = 649.5
#: Output cadence. Mercury's year is 88 days; 2-day samples resolve every
#: periodic term in its elements (the shortest that matter are ~44 days).
SAMPLE_DAYS = 2.0
#: Coarser cadence for the other planets (their shortest periods are months).
SAMPLE_DAYS_OUTER = 20.0

# --------------------------------------------------------------------------
# The bodies
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Body:
    key: str
    name: str
    #: NAIF ID of the DE440 segment (barycentric) holding this body.
    naif: int
    #: The GM constant's name in DE440's header.
    gm_key: str
    note: str = ""


#: The Sun, the planets as system barycentres, and Pluto. Earth and Moon travel
#: as one body at their barycentre: seen from Mercury, the pair's quadrupole is
#: (M_moon/M_earth) x (384,400 km / 0.6 au)^2 ~ 2e-7 of its pull, far below
#: anything measured here.
BODIES: list[Body] = [
    Body("sun", "Sun", 10, "GMS"),
    Body("mercury", "Mercury", 1, "GM1"),
    Body("venus", "Venus", 2, "GM2"),
    Body("emb", "Earth-Moon", 3, "GMB", "Earth-Moon barycentre"),
    Body("mars", "Mars", 4, "GM4", "system barycentre"),
    Body("jupiter", "Jupiter", 5, "GM5", "system barycentre"),
    Body("saturn", "Saturn", 6, "GM6", "system barycentre"),
    Body("uranus", "Uranus", 7, "GM7", "system barycentre"),
    Body("neptune", "Neptune", 8, "GM8", "system barycentre"),
    Body("pluto", "Pluto", 9, "GM9", "Pluto-Charon barycentre"),
]
BODY = {b.key: b for b in BODIES}
ALL_KEYS = [b.key for b in BODIES]
PLANET_KEYS = ALL_KEYS[1:]
PERTURBER_KEYS = ALL_KEYS[2:]  # everything that tugs on Mercury

# --------------------------------------------------------------------------
# The N-body experiments
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class RunSpec:
    """One numerical integration: which bodies, which forces, which integrator."""

    name: str
    bodies: tuple[str, ...]
    #: REBOUNDx effects: "gr" (1PN from the Sun only, Anderson et al. 1975),
    #: "gr_full" (1PN for every pair, Einstein-Infeld-Hoffmann), "gr_potential"
    #: (Nobili & Roxburgh's effective potential), "j2" (the Sun's oblateness),
    #: "lt" (the Sun's spin: Lense-Thirring).
    effects: tuple[str, ...] = ()
    integrator: str = "ias15"
    #: WHFast time step in days (ignored by IAS15).
    dt_days: float | None = None
    #: Multiply 1/c^2 by this (the exaggerated rosette uses 10^6).
    gr_boost: float = 1.0
    back_yr: float = SPAN_BACK_YR
    forward_yr: float = SPAN_FORWARD_YR
    description: str = ""
    #: Output cadence (days) for Mercury, and every how many samples the other
    #: planets are recorded.
    sample_days: float = SAMPLE_DAYS
    outer_every: int = int(round(SAMPLE_DAYS_OUTER / SAMPLE_DAYS))
    #: Rescale bodies' masses: ((key, factor), ...). Used to show which part of a
    #: planet's effect is first order in its mass.
    mass_scale: tuple = ()


TWO = ("sun", "mercury")
ALL = tuple(ALL_KEYS)

RUNS: list[RunSpec] = [
    # The numerical floor: a lone planet must not precess at all.
    RunSpec("two_body_newton", TWO, description="Sun + Mercury, Newton: must not precess"),
    # General relativity alone, three ways.
    RunSpec("two_body_gr", TWO, ("gr",), description="Sun + Mercury, 1PN (Anderson et al. 1975)"),
    RunSpec("two_body_gr_full", TWO, ("gr_full",), description="Sun + Mercury, 1PN EIH"),
    RunSpec("two_body_gr_potential", TWO, ("gr_potential",),
            description="Sun + Mercury, Nobili-Roxburgh potential"),
    # The Sun's shape and spin.
    RunSpec("two_body_j2", TWO, ("j2",), description="Sun + Mercury, solar oblateness J2"),
    RunSpec("two_body_lt", TWO, ("lt",), description="Sun + Mercury, solar Lense-Thirring"),
    # Each planet on its own.
    *[RunSpec(f"pair_{k}", ("sun", "mercury", k),
              description=f"Sun + Mercury + {BODY[k].name}, Newton") for k in PERTURBER_KEYS],
    # All of them together.
    RunSpec("all_newton", ALL, description="all planets, Newton"),
    RunSpec("all_gr", ALL, ("gr",), description="all planets, 1PN from the Sun"),
    RunSpec("all_gr_full", ALL, ("gr_full",), description="all planets, full 1PN (EIH)"),
    RunSpec("full_model", ALL, ("gr_full", "j2", "lt"),
            description="all planets, EIH + solar J2 + Lense-Thirring (closest to DE440)"),
]
RUN = {r.name: r for r in RUNS}

#: Leave one out: all planets but one. Park et al. (2017, Table 3) define each
#: planet's share this way (its mass set to zero in the full system), which keeps
#: the other planets' orbits realistic; the pair runs above do not (Jupiter
#: without Saturn follows a different average orbit).
LOO_RUNS: list[RunSpec] = [
    RunSpec(f"all_but_{k}", tuple(b for b in ALL if b != k),
            description=f"all planets except {BODY[k].name}, Newton")
    for k in ("venus", "emb", "mars", "jupiter", "saturn", "uranus", "neptune")
]
RUN.update({r.name: r for r in LOO_RUNS})

#: A symplectic integrator with too long a step fakes precession: WHFast runs of
#: the Newtonian Solar System at several steps, compared with IAS15. Shorter
#: spans are enough to measure the spurious rate.
WHFAST_STEPS_DAYS = (0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0)
WHFAST_SPAN_YR = 200.0
WHFAST_RUNS: list[RunSpec] = [
    RunSpec(f"whfast_dt{dt:g}", ALL, (), "whfast", dt, back_yr=0.0, forward_yr=WHFAST_SPAN_YR,
            description=f"all planets, Newton, WHFast dt = {dt:g} d")
    for dt in WHFAST_STEPS_DAYS
] + [RunSpec("ias15_200yr", ALL, (), back_yr=0.0, forward_yr=WHFAST_SPAN_YR,
             description="all planets, Newton, IAS15 (reference for WHFast)")]
RUN.update({r.name: r for r in WHFAST_RUNS})

#: Six millennia: our integration against JPL's DE441 (via Horizons) from
#: 1000 BC to AD 5000, far beyond DE440's span. The step divides 25 Julian
#: years exactly, so every Horizons epoch (one per 25 years) is a sample.
MILLENNIA_YR = 3000.0
MILLENNIA_SAMPLE_DAYS = 25.0 * YEAR_DAYS / 2300.0
MILLENNIA_RUNS: list[RunSpec] = [
    RunSpec(f"millennia_{tag}", ALL, eff, back_yr=MILLENNIA_YR, forward_yr=MILLENNIA_YR,
            sample_days=MILLENNIA_SAMPLE_DAYS, outer_every=10, description=desc)
    for tag, eff, desc in (("newton", (), "all planets, Newton, 1000 BC - AD 5000"),
                           ("full", ("gr_full", "j2", "lt"),
                            "all planets, EIH + solar J2 + Lense-Thirring, 1000 BC - AD 5000"))
]
RUN.update({r.name: r for r in MILLENNIA_RUNS})

#: Is Jupiter's share first order in its mass? The same pair run with Jupiter
#: 10 and 100 times lighter, scaled back up (Gauss's averaging is exactly
#: first order, the N-body run is not).
MASS_SCALE_RUNS: list[RunSpec] = [
    RunSpec(f"pair_jupiter_x{f:g}", ("sun", "mercury", "jupiter"), mass_scale=(("jupiter", f),),
            description=f"Sun + Mercury + Jupiter with Jupiter's mass x {f:g}")
    for f in (0.1, 0.01)
]
RUN.update({r.name: r for r in MASS_SCALE_RUNS})

#: Reversibility test: integrate forward this far, come back to J2000, compare.
REVERSIBILITY_YR = 649.5

# --------------------------------------------------------------------------
# The exaggerated animation
# --------------------------------------------------------------------------

#: General relativity turned up a million times: 1/c^2 x 10^6, i.e. light
#: slowed to 300 km/s. Mercury's perihelion then turns ~29 degrees per orbit.
ROSETTE_BOOST = 1.0e6
ROSETTE_ORBITS = 26

# --------------------------------------------------------------------------
# Published numbers, consulted only for comparison (each with its source)
# --------------------------------------------------------------------------

#: Physical constants used outside DE440 (IAU 2015 nominal / CODATA 2018).
GM_SUN_SI = 1.32712440041e20  # m^3 s^-2 (DE440: 132712440041.279 km^3/s^2)
GM_EARTH_SI = 3.986004418e14  # m^3 s^-2 (IERS)
C_SI = 299792458.0
AU_M = 149597870700.0
PC_M = 3.0856775814913673e16
R_SUN_M = 6.957e8  # IAU 2015 nominal solar radius
T_SUN_S = GM_SUN_SI / C_SI**3  # 4.925490947 microseconds
JULIAN_YEAR_S = 365.25 * 86400.0

PUBLISHED = {
    "park2017": {
        "total_arcsec_cy": 575.3100, "total_err": 0.0015,
        "beta_minus_1": -2.7e-5, "beta_minus_1_err": 3.9e-5,
        "j2_sun": 2.25e-7, "j2_sun_err": 0.09e-7,
        "source": "Park et al. 2017, AJ 153, 121 (MESSENGER ranging)",
    },
    "cassini_gamma": {"gamma_minus_1": 2.1e-5, "err": 2.3e-5,
                      "source": "Bertotti, Iess & Tortora 2003, Nature 425, 374"},
    #: GRAVITY's Schwarzschild-precession factor for S2 (GR = 1, Newton = 0).
    "gravity_s2": [
        {"label": "GRAVITY 2020", "f_sp": 1.10, "err": 0.19, "stars": "S2",
         "source": "GRAVITY Collaboration 2020, A&A 636, L5 (Table E.1)"},
        {"label": "GRAVITY 2022", "f_sp": 0.85, "err": 0.16, "stars": "S2",
         "source": "GRAVITY Collaboration 2022, A&A 657, L12 (Sect. 4)"},
        {"label": "GRAVITY 2022", "f_sp": 0.997, "err": 0.144, "stars": "S2, S29, S38, S55",
         "source": "GRAVITY Collaboration 2022, A&A 657, L12 (Table B.1)"},
        {"label": "GRAVITY 2024", "f_sp": 0.918, "err": 0.128, "stars": "S2",
         "source": "GRAVITY Collaboration 2024, arXiv:2409.12261v1 (Sect. 3.2; submitted version)"},
        {"label": "GRAVITY 2024", "f_sp": 1.135, "err": 0.110, "stars": "S2, S29, S38, S55",
         "source": "GRAVITY Collaboration 2024, arXiv:2409.12261v1 (Table B.1; submitted version)"},
    ],
    #: S2's orbit and Sgr A* (GRAVITY 2020, Table E.1; osculating at 2010.35).
    "s2_orbit_gravity2020": {
        "a_mas": 125.058, "e": 0.884649, "inc_deg": 134.567, "omega_deg": 66.263,
        "Omega_deg": 228.171, "period_yr": 16.0455, "t_peri_yr": 2018.37900,
        "mass_1e6": 4.261, "r0_pc": 8246.7,
        "source": "GRAVITY Collaboration 2020, A&A 636, L5, Table E.1"},
    #: PSR B1913+16 (Weisberg & Huang 2016, ApJ 829, 55, Table 2 and Sect. 4).
    #: gamma is printed "0.004307 ms" in the paper: it is 4.307 ms (the TEMPO
    #: parameter file shipped with the arXiv source has GAMMA 0.004306748 s).
    "b1913": {
        "pb_days": 0.322997448918, "e": 0.6171340, "x_s": 2.341776,
        "omdot_deg_yr": 4.226585, "omdot_err": 0.000004,
        "gamma_s": 4.307e-3, "gamma_err": 0.004e-3,
        "pbdot_intrinsic": -2.398e-12, "pbdot_intrinsic_err": 0.004e-12,
        "m_p": 1.438, "m_c": 1.390, "pbdot_ratio": 0.9983, "pbdot_ratio_err": 0.0016,
        "source": "Weisberg & Huang 2016, ApJ 829, 55"},
    #: The double pulsar (Kramer et al. 2021, PRX 11, 041050, Table IV/V, Eqs. 36-47);
    #: x_B and R from Kramer et al. 2006, Science 314, 97 (Table 1).
    "j0737": {
        "pb_days": 0.1022515592973, "e": 0.087777023, "x_a_s": 1.415028603,
        "omdot_deg_yr": 16.899323, "omdot_err": 0.000013,
        "gamma_s": 0.384045e-3, "gamma_err": 0.000094e-3,
        "s": 0.999936, "s_err_plus": 0.000009, "s_err_minus": 0.000010,
        "mass_ratio": 1.0714, "mass_ratio_err": 0.0011,
        "pbdot_gw": -1.247782e-12, "pbdot_gw_err": 0.000079e-12,
        "m_a": 1.338185, "m_b": 1.248868,
        "k_test": 1.000015, "k_test_err": 0.000026,
        "source": "Kramer et al. 2021, PRX 11, 041050; R from Kramer et al. 2006"},
    #: LAGEOS II perigee (Lucchesi & Peron 2010, PRL 105, 231103): epsilon
    #: scales the whole relativistic perigee rate (Schwarzschild 3351.95 mas/yr,
    #: Lense-Thirring -57.00, de Sitter +10.69). The authors add a preliminary
    #: ~2% systematic error (mostly the Earth's even zonal harmonics).
    "lageos2": {"a_km": 12162.0, "e": 0.014, "eps": 1.00028, "eps_stat": 0.00214,
                "eps_sys": 0.02, "schwarzschild_mas_yr": 3351.95,
                "source": "Lucchesi & Peron 2010, PRL 105, 231103"},
    #: (1566) Icarus: lambda = 0.95 +- 0.08 (formal); the authors' realistic error ~0.2.
    "icarus": {"a_au": 1.078, "e": 0.827, "ratio": 0.95, "err_formal": 0.08, "err": 0.2,
               "source": "Shapiro, Smith, Ash & Herrick 1971, AJ 76, 588"},
    #: DI Herculis: rates in degrees per orbital cycle.
    "di_her": {
        "p_days": 10.5501696, "e": 0.4895, "m1": 5.17, "m2": 4.524, "r1": 2.681, "r2": 2.478,
        "obliquity_1_deg": 75.0, "obliquity_2_deg": 80.0,
        "history": [
            {"label": "Guinan & Maloney 1985 (all minima 1900-84)", "year": 1985,
             "obs": 0.000188, "err": 0.000052},
            {"label": "Guinan et al. 1994 (IBVS 4101)", "year": 1994, "obs": 0.000300,
             "err": 0.000043},
            {"label": "Claret, Torres & Wolf 2010", "year": 2010, "obs": 0.00042, "err": 0.00012},
            {"label": "Claret et al. 2021 (TESS)", "year": 2021, "obs": 0.000435, "err": 0.000023},
            {"label": "Liang, Winn & Albrecht 2022 (TESS)", "year": 2022, "obs": 0.000392,
             "err_plus": 0.000108, "err_minus": 0.000078, "err": 0.000093},
        ],
        "gr_gm85": 0.000676, "classical_aligned_gm85": 0.000557, "classical_aligned_gm85_err": 0.000075,
        "theory_aligned_gm85": 0.001233, "theory_aligned_gm85_err": 0.000087,
        "gr_claret10": 0.000677, "tidal_claret10": 0.000339, "rotation_misaligned_claret10": -0.00056,
        "theory_misaligned_claret21": 0.00046, "theory_misaligned_claret21_err": 0.00009,
        "theory_misaligned_liang22": 0.000375, "theory_misaligned_liang22_err": 0.000150,
        "sources": ("Guinan & Maloney 1985, AJ 90, 1519; Guinan et al. 1994, IBVS 4101; "
                    "Albrecht et al. 2009, Nature 461, 373; Claret, Torres & Wolf 2010, A&A 515, A4; "
                    "Claret et al. 2021, A&A 654, A17; Liang, Winn & Albrecht 2022, ApJ 927, 114"),
    },
    #: EW Ori, added in paper II (Claret et al. 2021, Appendix A).
    "ew_ori": {"wdot_rel_mea": 0.00026, "err": 0.00002, "obs": 0.00033, "obs_err": 0.00002,
               "source": "Claret et al. 2021, A&A 654, A17, Appendix A"},
    #: Baroch et al.'s combined test: slope A = (2 + 2 gamma - beta)/3.
    "baroch2021_A": {"A": 1.002, "err": 0.012, "A_paper2": 1.000, "err_paper2": 0.011,
                     "source": "Baroch et al. 2021, A&A 649, A64 (Sect. 6); paper II: Claret et al. 2021"},
}

#: Mercury's anomalous (non-Newtonian) perihelion advance through history, in
#: "/century. "primary": read in the original (most via the N7 project's
#: literature checks, 64-gravity-scorecard/literature/checks/beta_lunar.csv;
#: Clemence 1943 and Anderson et al. 1991 from the ADS page scans);
#: "secondary": known only from a later compilation, named in ``source``.
#: Errors are 1 sigma where the original says so; probable errors x 1.4826.
MERCURY_HISTORY = [
    {"year": 1859, "who": "Le Verrier", "value": 38.0, "err": None, "quality": "primary",
     "source": "Le Verrier 1859, Comptes rendus 49, 379 (p. 380): '38 secondes'",
     "note": "transits of Mercury against his planetary theory; no error given"},
    {"year": 1882, "who": "Newcomb", "value": 42.95, "err": None, "quality": "primary",
     "source": "Newcomb 1882, Astron. Papers Am. Ephem. 1, 363 (p. 473: p = +42.95)",
     "note": "transits 1677-1881"},
    {"year": 1895, "who": "Newcomb", "value": 41.24, "err": 2.09, "quality": "primary",
     "source": "Newcomb 1895, The Elements of the Four Inner Planets, p. 109: e x excess = 8.48 +- 0.43",
     "note": "divided by e = 0.2056; 'mean errors' (1 sigma). Often misquoted as 43\""},
    {"year": 1943, "who": "Clemence", "value": 42.84, "err": 0.61, "quality": "primary",
     "source": "Clemence 1943, AJ 50, 126 (p. 127): O - C = -0.07 +- 0.41 (p.e.) after the relativity effect",
     "note": "relativity term 42.91 = Newcomb's empirical 43.37 (secondary) - 0.46 (p. 127); "
             "error excludes the planetary masses"},
    {"year": 1947, "who": "Clemence", "value": 42.56, "err": 0.94, "quality": "secondary",
     "source": "Clemence 1947, Rev. Mod. Phys. 19, 361, via Pireaux & Rozelot 2003 (Table 1 notes)",
     "note": "5599.74 +- 0.41 observed minus 5557.18 +- 0.85 computed; error type unverified"},
    {"year": 1958, "who": "Duncombe", "value": 43.10, "err": 0.44, "quality": "secondary",
     "source": "Duncombe 1958, via Pireaux & Rozelot 2003, Table 1", "note": "optical meridian + transits"},
    {"year": 1972, "who": "Shapiro et al. (radar)", "value": 1.005 * 42.98, "err": 0.02 * 42.98,
     "quality": "primary",
     "source": "Shapiro et al. 1972, PRL 28, 1594: (2+2g-b)/3 = 1.005 +- 0.007 (stat), ~0.02 with systematics",
     "note": "ratio to GR x 42.98; the authors' realistic error"},
    {"year": 1976, "who": "Shapiro et al. (radar)", "value": 43.11, "err": 0.21, "quality": "secondary",
     "source": "via Anderson et al. 1991 (Proc. ASA 9, 324) and Pireaux & Rozelot 2003",
     "note": "attributed to Shapiro, Counselman & King 1976, PRL 36, 555 (paywalled; the abstract gives only beta and gamma)"},
    {"year": 1987, "who": "Anderson et al. (radar + Mariner 10)", "value": 42.92, "err": 0.20,
     "quality": "secondary", "source": "Anderson et al. 1987, Icarus 71, 337, via Anderson et al. 1991", "note": ""},
    {"year": 1991, "who": "Anderson et al. (radar 1966-88)", "value": 42.94, "err": 0.20, "quality": "primary",
     "source": "Anderson et al. 1991, Proc. ASA 9, 324: 'an excess precession of 42.94 arcsec'",
     "note": "the authors could not improve on the 0.20 error of 1987"},
    {"year": 1992, "who": "Anderson et al.", "value": 43.13, "err": 0.14, "quality": "secondary",
     "source": "Anderson et al. 1992, via Pireaux & Rozelot 2003, Table 1", "note": ""},
    {"year": 2001, "who": "Pitjeva (EPM2000)", "value": 43.0115, "err": 0.0085, "quality": "secondary",
     "source": "Pitjeva 2001, via Pireaux & Rozelot 2003, Table 1", "note": "planetary ephemeris fit"},
    {"year": 2017, "who": "Park et al. (MESSENGER)", "value": 42.9799, "err": 0.0009, "quality": "primary",
     "source": "Park et al. 2017, AJ 153, 121, Table 3 (gravitoelectric term)",
     "note": "with Cassini's gamma and their fitted beta and J2"},
]

#: Park et al. 2017, Table 3: contributions to the in-plane precession rate in
#: Mercury's mean orbital plane, "/Julian century. Planet rows are each planet's
#: direct effect (its GM set to zero); the 'interaction' rows are indirect effects
#: through the other planets' orbits.
PARK2017_TABLE3 = {
    "planets": {"venus": 277.4176, "emb": 90.8881, "mars": 2.4814, "jupiter": 153.9899,
                "saturn": 7.3227, "uranus": 0.1425, "neptune": 0.0424},
    "interactions": {"mercury": 0.0050, "mercury+venus": -0.0053, "venus+emb": -0.0209,
                     "venus+jupiter": -0.0012, "emb+mars": -0.0016, "mars+jupiter": 0.0002,
                     "jupiter+saturn": 0.0411, "saturn+uranus": 0.0004},
    "asteroids": 0.0012, "solar_oblateness": 0.0286, "solar_oblateness_err": 0.0011,
    "gravitoelectric": 42.9799, "gravitoelectric_err": 0.0009,
    "lense_thirring": -0.0020, "lense_thirring_err": 0.0002, "total": 575.3100, "total_err": 0.0015,
    "Q_arcsec_per_cy2": -0.04478132,
    "source": "Park et al. 2017, AJ 153, 121, Table 3 and Sect. 4",
}
