"""Published measurements and predictions that this project compares against.

Every number was checked against the paper itself (text extracts in
``data/literature/extracted/``); uncertainties are 1 sigma unless noted.
"""

from __future__ import annotations

import datetime as dt

# ---------------------------------------------------------------------------
# Gravity Probe B
# ---------------------------------------------------------------------------

#: Everitt et al. (2011), PRL 106, 221101, Table II; identical in Everitt et al.
#: (2015), CQG 32, 224001. mas/yr. NS = north-south (geodetic), WE = west-east
#: (frame dragging). ê_WE = ê_GS x z (local west at the guide star), ê_NS = ê_WE x ê_GS
#: (local north); negative NS = southward, negative WE = eastward.
GPB_RESULTS = {
    "Gyro 1": (-6588.6, 31.7, -41.3, 24.6),
    "Gyro 2": (-6707.0, 64.1, -16.1, 29.7),
    "Gyro 3": (-6610.5, 43.2, -25.0, 12.1),
    "Gyro 4": (-6588.7, 33.2, -49.3, 11.4),
    "Joint": (-6601.8, 18.3, -37.2, 7.2),
}
GPB_GR_PREDICTION = {"ns": -6606.1, "we": -39.2}   # Earth's geodetic and frame dragging only
#: Terms removed before comparing (PRL 2011): guide-star proper motion and the
#: solar geodetic effect (19.2 mas/yr about the ecliptic pole; Conklin et al. 2015).
GPB_REMOVED = {"proper_motion": {"ns": 27.3, "we": -20.0},
               "solar_geodetic": {"ns": 7.3, "we": -16.2}}
#: Joint error budget (NS / WE): statistical 16.8 / 5.9; parameter sensitivity 7.1 / 4.0.
GPB_ERROR_BUDGET = {"statistical": (16.8, 5.9), "parameter_sensitivity": (7.1, 4.0),
                    "guide_star_motion": (0.1, 0.1), "solar_geodetic": (0.3, 0.6),
                    "telescope_readout": (0.5, 0.5)}
#: The Earth's flattening reduces the geodetic rate by ~7 mas/yr for GP-B's polar
#: orbit (Breakwell 1988; "a 7 marc-s/yr correction to Omega_g", Everitt et al. 2015).
GPB_J2_CORRECTION = -7.0

#: Orbit (Stanford GP-B Post-Flight Analysis Final Report 2007, "Quick facts"):
#: apogee/perigee altitude 659.1/639.5 km.
GPB_ORBIT = {"a_m": 7_027_400.0, "e": 0.0014, "i_deg": 90.007, "argp_deg": 71.3,
             "node_deg": 163.26, "period_min": 97.6}
GPB_SCIENCE = (dt.date(2004, 8, 28), dt.date(2005, 8, 14))

#: Guide star IM Pegasi: VLBI astrometry (Bartel et al. 2015, CQG 32, 224021, Table 3;
#: Ratner et al. 2012, ApJS 201, 5), ICRF2, epoch 2005.08.
IM_PEG = {"ra_deg": (22 + 53 / 60 + 2.258612 / 3600) * 15,
          "dec_deg": 16 + 50 / 60 + 28.16005 / 3600,
          "epoch_mjd": 53_399.0,                       # 2005.08
          "pmra_mas_yr": -20.833, "pmdec_mas_yr": -27.267,   # mu_alpha* , mu_delta
          "parallax_mas": 10.370}

# ---------------------------------------------------------------------------
# Frame dragging with laser-ranged satellites
# ---------------------------------------------------------------------------

#: mu = measured / GR. (stat, syst) 1 sigma unless noted.
FRAME_DRAGGING_SLR = [
    {"ref": "Ciufolini & Pavlis 2004", "journal": "Nature 431, 958", "year": 2004,
     "sats": "LAGEOS + LAGEOS 2", "mu": 0.99, "stat": 0.05, "syst": 0.0,
     "note": "99 +- 5%, authors allow +-10% in total; 11 yr, EIGEN-GRACE02S", "total": 0.10},
    {"ref": "Ciufolini et al. 2016", "journal": "EPJC 76, 120", "year": 2016,
     "sats": "LARES + LAGEOS + LAGEOS 2", "mu": 0.994, "stat": 0.002, "syst": 0.05,
     "note": "Feb 2012 - Sep 2015, GGM05S"},
    {"ref": "Ciufolini et al. 2019", "journal": "EPJC 79, 872", "year": 2019,
     "sats": "LARES + LAGEOS + LAGEOS 2", "mu": 0.9910, "stat": 0.0006, "syst": 0.02,
     "note": "7 yr LARES, 26 yr LAGEOS; systematic up to 0.04"},
    {"ref": "Lucchesi et al. 2020", "journal": "Universe 6, 139", "year": 2020,
     "sats": "LAGEOS + LAGEOS II + LARES", "mu": 1.0015, "stat": 0.0074 / 1.96, "syst": 0.016,
     "note": "mu-1 = (1.5 +- 7.4) e-3 (95%) +- 16e-3; 2359 days from April 2012"},
    {"ref": "Ciufolini et al. 2026", "journal": "Nature 655, 332", "year": 2026,
     "sats": "LARES 2 + LAGEOS", "mu": 1.001, "stat": 0.001, "syst": 0.002,
     "note": "GEODYN II, GGM05S + GRACE time variation, ~3 yr (one 1050-day nodal period)"},
]
#: Lense-Thirring node rates (mas/yr) quoted in the papers.
LT_NODE_RATE = {"lageos1": 30.67, "lares2": 30.71, "lageos2": 31.50}
#: de Sitter (geodetic) precession of every orbit around the ecliptic pole.
DE_SITTER_RATE = 19.2           # mas/yr; on the equatorial node ~17.6 (x cos 23.4 deg)

#: Published orbits (Ciufolini et al. 2026, Extended Data Table 1; LARES 2 also
#: Ciufolini et al. 2023, EPJ Plus 138, 1054: a = 12,266.1359 km, i = 70.1615 deg).
SLR_ORBITS = {"lageos1": {"a_km": 12_270.020, "e": 0.004, "i_deg": 109.85},
              "lares2": {"a_km": 12_266.136, "e": 0.00027, "i_deg": 70.1615}}

#: GP-B's frame-dragging and geodetic results for the frame-dragging ladder
GPB_FRAME_DRAGGING_RATIO = (-37.2 / -39.2, 7.2 / 39.2)
GPB_GEODETIC_RATIO = (-6601.8 / -6606.1, 18.3 / 6606.1)

# ---------------------------------------------------------------------------
# Spin (geodetic) precession of pulsars, deg/yr
# ---------------------------------------------------------------------------

#: Masses (Msun) of the precessing star and its companion, orbital period (d), e.
PULSAR_SYSTEMS = {
    "J0737-3039B": {"pb_d": 0.1022515592973, "e": 0.087777023, "mp": 1.248868, "mc": 1.338185,
                    "ref": "Kramer et al. 2021, PRX 11, 041050"},
    "B1913+16": {"pb_d": 0.322997448918, "e": 0.6171340, "mp": 1.438, "mc": 1.390,
                 "ref": "Weisberg & Huang 2016, ApJ 829, 55"},
    "J1906+0746": {"pb_d": 0.16599304683, "e": 0.0853028, "mp": 1.291, "mc": 1.322,
                   "ref": "van Leeuwen et al. 2015, ApJ 798, 118"},
    "B1534+12": {"pb_d": 0.420737298881, "e": 0.27367740, "mp": 1.3330, "mc": 1.3455,
                 "ref": "Fonseca et al. 2014, ApJ 787, 82"},
    "J1141-6545": {"pb_d": 0.19765096149, "e": 0.171876, "mp": 1.27, "mc": 1.02,
                   "ref": "Venkatraman Krishnan et al. 2020, Science 367, 577"},
}
#: Measured spin-precession rates (deg/yr): value, +err, -err, published GR value.
PULSAR_PRECESSION = [
    {"system": "J0737-3039B", "ref": "Breton et al. 2008", "journal": "Science 321, 104",
     "value": 4.77, "plus": 0.66, "minus": 0.65, "gr": 5.0734},
    {"system": "J0737-3039B", "ref": "Lower et al. 2024", "journal": "A&A 682, A26",
     "value": 5.16, "plus": 0.32, "minus": 0.34, "gr": 5.074005},
    {"system": "J1906+0746", "ref": "Desvignes et al. 2019", "journal": "Science 365, 1013",
     "value": 2.17, "plus": 0.11, "minus": 0.11, "gr": 2.234},
    {"system": "B1534+12", "ref": "Stairs et al. 2004", "journal": "PRL 93, 141101",
     "value": 0.44, "plus": 0.48, "minus": 0.16, "gr": 0.51},
    {"system": "B1534+12", "ref": "Fonseca et al. 2014", "journal": "ApJ 787, 82",
     "value": 0.59, "plus": 0.12, "minus": 0.08, "gr": 0.51},
]
#: No measured rate: B1913+16's profile changes (Kramer 1998; Weisberg & Taylor 2002)
#: and J1141-6545's (Hotan et al. 2005; Manchester et al. 2010) are fitted assuming
#: GR's rate. J1141-6545's 2020 result is a different effect: its orbit is tilted
#: by the spin of its fast-rotating white dwarf (Lense-Thirring plus the white
#: dwarf's spin-induced quadrupole), from x-dot = (1.7 +- 0.3) e-13.
PULSAR_NOT_MEASURED = {"B1913+16": "Kramer 1998, ApJ 509, 856; Weisberg & Taylor 2002, ApJ 576, 942",
                       "J1141-6545": "Hotan et al. 2005; Manchester et al. 2010; "
                                     "Venkatraman Krishnan et al. 2019, ApJL 873, L15"}

# ---------------------------------------------------------------------------
# Lunar laser ranging: the de Sitter precession of the Moon's orbit (19.2 mas/yr)
# ---------------------------------------------------------------------------

#: Fractional precision of the geodetic-precession test (deviation 0 assumed unless given).
LLR_DE_SITTER = [
    {"ref": "Bertotti, Ciufolini & Bender 1987", "journal": "PRL 58, 1062", "year": 1987,
     "dev": 0.0, "err": 0.10},
    {"ref": "Shapiro et al. 1988", "journal": "PRL 61, 2643", "year": 1988, "dev": 0.0, "err": 0.02},
    {"ref": "Williams, Newhall & Dickey 1996", "journal": "PRD 53, 6730", "year": 1996,
     "dev": 0.0, "err": 0.007},
    {"ref": "Williams, Turyshev & Boggs 2004", "journal": "PRL 93, 261101", "year": 2004,
     "dev": -0.0019, "err": 0.0064},
    {"ref": "Hofmann & Müller 2018", "journal": "CQG 35, 035015", "year": 2018,
     "dev": 0.0, "err": 0.0009},
]
