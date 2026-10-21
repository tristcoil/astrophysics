"""Published numbers: the census of known thermonuclear runaways, and the
values we grade our own analysis against.

Everything here was copied by hand from the papers (arXiv versions read on
4 Oct 2026) and is used for two things only: the radial velocities, which no
public catalogue holds for these faint, peculiar stars, and comparisons.
Positions, parallaxes and proper motions always come from Gaia DR3 itself.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Runaway:
    name: str
    source_id: int            # Gaia DR3
    cls: str                  # "D6", "LP40-365", "He-star donor", "other"
    rv: float                 # heliocentric radial velocity, km/s
    rv_err: float
    rv_ref: str
    discovery: str
    #: El-Badry et al. (2023) Table 3, fiducial (flat-M_G) prior: median with
    #: (minus, plus) 68% half-widths. None where they list no value.
    eb23: dict = field(default_factory=dict)
    note: str = ""


def _v(med, lo, hi):
    return (med, lo, hi)


#: Every object credibly proposed as a runaway from a thermonuclear supernova,
#: as of October 2026: El-Badry et al. (2023, Table 3: 16 objects), plus
#: J1022-3414 (Bhat et al. 2026). J1637+3631's velocity is Hollands et al.
#: (2025)'s spectrum, which supersedes the SDSS template estimate (300 +- 50)
#: in El-Badry et al.'s table.
CENSUS: list[Runaway] = [
    # --- "D6 stars": surviving donors of double-degenerate double detonations
    Runaway("D6-1", 5805243926609660032, "D6", 1200, 40, "Shen et al. 2018", "Shen et al. 2018",
            eb23=dict(d=_v(1.91, 0.22, 0.29), vtot=_v(2045, 187, 251), vej=_v(2254, 185, 248),
                      z=_v(-0.57, 0.09, 0.07), tdisk=_v(0.63, 0.04, 0.04))),
    Runaway("D6-2", 1798008584396457088, "D6", 80, 10, "Chandra et al. 2022 (via El-Badry et al. 2023)",
            "Shen et al. 2018",
            eb23=dict(d=_v(0.84, 0.04, 0.05), vtot=_v(1151, 51, 59), vej=_v(1051, 54, 62),
                      z=_v(-0.27, 0.02, 0.01), tdisk=_v(-0.72, 0.01, 0.01)),
            note="LP 398-9; Shen et al. (2018) measured RV = 20 +- 60 km/s"),
    Runaway("D6-3", 2156908318076164224, "D6", -20, 80, "Shen et al. 2018", "Shen et al. 2018",
            eb23=dict(d=_v(2.26, 0.35, 0.34), vtot=_v(2248, 353, 340), vej=_v(2393, 391, 377),
                      z=_v(0.93, 0.14, 0.14), tdisk=_v(2.24, 0.17, 0.20))),
    Runaway("J1235-3752", 6156470924553703552, "D6", -1694, 10, "El-Badry et al. 2023",
            "El-Badry et al. 2023",
            eb23=dict(d=_v(4.07, 1.19, 1.01), vtot=_v(2670, 350, 339), vej=_v(2471, 345, 351),
                      z=_v(1.73, 0.50, 0.42), tdisk=_v(1.87, 0.25, 0.86))),
    Runaway("J0927-6335", 5250394728194220800, "D6", -2285, 20, "El-Badry et al. 2023",
            "El-Badry et al. 2023",
            eb23=dict(d=_v(4.53, 1.71, 2.17), vtot=_v(2753, 148, 271), vej=_v(2519, 147, 271),
                      z=_v(-0.70, 0.35, 0.27), tdisk=_v(-1.32, 0.42, 0.42))),
    Runaway("J0546+0836", 3335306915849417984, "D6", 1200, 20, "El-Badry et al. 2023",
            "El-Badry et al. 2023",
            eb23=dict(d=_v(4.02, 1.67, 2.25), vtot=_v(1699, 390, 670), vej=_v(1864, 416, 682),
                      z=_v(-0.69, 0.40, 0.29), tdisk=_v(0.61, 0.08, 0.05))),
    Runaway("J1332-3541", 6164642052589392512, "D6", 1090, 50, "El-Badry et al. 2023",
            "El-Badry et al. 2023",
            eb23=dict(d=_v(1.63, 0.68, 1.19), vtot=_v(1464, 331, 740), vej=_v(1619, 320, 707),
                      z=_v(0.74, 0.30, 0.53), tdisk=_v(0.64, 0.14, 0.13))),
    Runaway("J1637+3631", 1327920737357113088, "D6", 384, 6, "Hollands et al. 2025",
            "Raddi et al. 2019 (classified D6 by Hollands et al. 2025)",
            eb23=dict(d=_v(3.28, 1.74, 3.18), vtot=_v(1202, 463, 962), vej=_v(1030, 467, 1205),
                      z=_v(2.20, 1.16, 2.12), tdisk=_v(4.26, 1.33, 1.03)),
            note="El-Badry et al. listed it as 'other' with RV 300 +- 50 from the SDSS spectrum; "
                 "Hollands et al. (2025): ejection 1870 (+360 -300) km/s, flight time 4.5 Myr"),
    # --- "LP 40-365 stars": partly burnt accretors, probably from Type Iax events
    Runaway("LP 40-365", 1711956376295435520, "LP40-365", 498, 5, "Raddi et al. 2019",
            "Vennes et al. 2017",
            eb23=dict(d=_v(0.61, 0.01, 0.01), vtot=_v(837, 5, 5), vej=_v(607, 5, 5),
                      z=_v(0.43, 0.01, 0.01), tdisk=_v(4.87, 0.33, 0.39)), note="GD 492"),
    Runaway("J1603-6613", 5822236741381879040, "LP40-365", -480, 5, "Raddi et al. 2019",
            "Raddi et al. 2019",
            eb23=dict(d=_v(2.07, 0.32, 0.49), vtot=_v(833, 38, 62), vej=_v(606, 42, 73),
                      z=_v(-0.35, 0.09, 0.06), tdisk=_v(1.58, 0.10, 0.12))),
    Runaway("J0905+2510", 688380457508503040, "LP40-365", 300, 50, "Raddi et al. 2019",
            "Raddi et al. 2019",
            eb23=dict(d=_v(4.49, 2.76, 7.59), vtot=_v(519, 271, 1122), vej=_v(737, 341, 1123),
                      z=_v(2.90, 1.77, 4.86), tdisk=_v(-13.58, 15.41, 32.35))),
    Runaway("J1825-3757", 6727110900983876096, "LP40-365", -47, 5, "Raddi et al. 2019",
            "Raddi et al. 2019",
            eb23=dict(d=_v(0.95, 0.03, 0.03), vtot=_v(429, 17, 18), vej=_v(662, 17, 18),
                      z=_v(-0.17, 0.01, 0.01), tdisk=_v(1.60, 0.02, 0.02))),
    Runaway("J1311-1846", 3507697866498687232, "LP40-365", 55, 10, "El-Badry et al. 2023",
            "El-Badry et al. 2023",
            eb23=dict(d=_v(1.88, 0.51, 1.09), vtot=_v(952, 201, 430), vej=_v(827, 217, 461),
                      z=_v(1.32, 0.35, 0.76), tdisk=_v(3.88, 0.20, 0.20))),
    Runaway("J1109+0001", 3804182280735442560, "LP40-365", 100, 10, "El-Badry et al. 2023",
            "El-Badry et al. 2023",
            eb23=dict(d=_v(2.90, 1.24, 2.31), vtot=_v(1378, 510, 957), vej=_v(1171, 510, 993),
                      z=_v(2.35, 0.99, 1.86), tdisk=_v(3.92, 0.38, 0.26))),
    Runaway("J1022-3414", 5446737753669901568, "LP40-365", -390, 20, "Bhat et al. 2026",
            "Bhat et al. 2026",
            note="ejected ~2.8 Myr ago at > 600 km/s (Bhat et al. 2026)"),
    # --- a helium-star donor, and one oxygen-atmosphere white dwarf
    Runaway("US 708", 815106177700219392, "He-star donor", 917, 7, "Geier et al. 2015",
            "Hirsch et al. 2005",
            eb23=dict(d=_v(8.38, 1.01, 1.00), vtot=_v(994, 10, 10), vej=_v(833, 16, 17),
                      z=_v(6.17, 0.74, 0.74), tdisk=_v(11.39, 1.70, 1.82)),
            note="hot subdwarf; distance from its spectrum: 8.5 +- 1.0 kpc (Geier et al. 2015)"),
    Runaway("J1240+6710", 1682129610835350400, "other", -177, 10, "Gänsicke et al. 2020",
            "Kepler et al. 2016",
            eb23=dict(d=_v(0.42, 0.02, 0.02), vtot=_v(239, 18, 20), vej=_v(446, 18, 21),
                      z=_v(0.35, 0.02, 0.02)),
            note="oxygen-atmosphere white dwarf, a partly burnt remnant"),
]

CLASS_LABELS = {
    "D6": "D⁶ stars (surviving donors)",
    "LP40-365": "LP 40-365 stars (Type Iax remnants?)",
    "He-star donor": "US 708 (helium-star donor)",
    "other": "J1240+6710 (oxygen white dwarf)",
}

#: US 708's spectroscopic distance (Geier et al. 2015), used as the prior for
#: its distance in place of the flat-magnitude prior, as El-Badry et al. did.
US708_DISTANCE_PRIOR_KPC = (8.5, 1.0)

#: El-Badry et al. (2023), Table 1: the 25 sources their query returned,
#: in their order, with their follow-up verdicts.
ELBADRY_TABLE1 = [
    (6156470924553703552, "hot D6 star (J1235-3752)"),
    (2156908318076164224, "D6 star (D6-3)"),
    (3335306915849417984, "hot D6 star (J0546+0836)"),
    (5250394728194220800, "hot D6 star (J0927-6335)"),
    (6164642052589392512, "hot D6 star (J1332-3541)"),
    (3542263595793124480, "bright neighbour"),
    (3804182280735442560, "LP 40-365 star (J1109+0001)"),
    (5703888058542880896, "DA white dwarf"),
    (5517276097516408576, "bright neighbour"),
    (4771417432717575680, "DA white dwarf"),
    (4546525523392712064, "DA white dwarf"),
    (6853349473073333632, "DA white dwarf"),
    (2393804867149529856, "DA white dwarf"),
    (3507697866498687232, "LP 40-365 star (J1311-1846)"),
    (5998866829060560768, "main-sequence star"),
    (5183592902806605824, "not observed"),
    (4129413800145771776, "DA white dwarf"),
    (1512757030058082304, "DA white dwarf"),
    (1271663056690275712, "DA white dwarf"),
    (6536783647184624640, "DA white dwarf"),
    (4096984529352659584, "main-sequence star"),
    (6592388973858801920, "DA white dwarf"),
    (6688913592127235584, "DA white dwarf"),
    (4096331041464263296, "main-sequence star"),
    (6703717691563155968, "DA white dwarf"),
]

#: Shen et al. (2018), Table 1: the seven Gaia DR2 candidates they followed
#: up. O1-O4 turned out to be ordinary hydrogen-rich stars, "whose true
#: parallaxes are likely much larger than measured". (name, DR2 source_id,
#: DR2 parallax, error, mas).
SHEN18_TABLE1 = [
    ("D6-1", 5805243926609660032, 0.471, 0.102),
    ("D6-2", 1798008584396457088, 1.052, 0.109),
    ("D6-3", 2156908318076164224, 0.427, 0.126),
    ("O1", 2050179518946705152, 0.237, 0.317),
    ("O2", 4396109004117478656, 2.186, 0.225),
    ("O3", 5884527618445501056, 0.360, 0.571),
    ("O4", 1820931585123817728, 0.574, 0.076),
]

#: D6-2 and its remnant (Shen et al. 2018, section 4.4): 9e4 yr ago the star
#: was 0.9 deg from the centre of G70.0-21.5 (Fesen et al. 2015), whose
#: distance from its shock velocities is 1-2 kpc.
D62_REMNANT = dict(name="G070.0-21.5", lookback_yr=9e4, offset_deg=0.9, distance_kpc=(1.0, 2.0))

#: Type Ia remnants searched by Shen et al. (2018) for survivors.
TYPE_IA_REMNANTS = ["G004.5+06.8", "G120.1+01.4", "G315.4-02.3", "G327.6+14.6"]

#: Published comparison numbers.
REFERENCE = {
    # Shen et al. 2018, section 2.4: runaway speeds for a 1.0 Msun accretor
    "shen18_vrun_0.6+1.0": 1800.0,
    "shen18_vrun_1.0+1.0": 2200.0,
    # El-Badry et al. 2023, section 7 (Macc <= 1.05 Msun, Bedard et al. radii)
    "eb23_vmax_two_1.05": 2235.0,
    "eb23_J0927_vej_2sigma_lower": 2315.0,
    # Birth-rate estimates, as fractions of the Type Ia rate
    "eb23_rate_fast": (0.03, 0.05),     # the five fast D6 stars
    "eb23_rate_slow": 1.0,              # D6-2 and J1332-like stars, ~100%
    "shen25_fraction": 0.02,            # Shen 2025: ~2% leave a survivor
    "igoshev23": "at least two orders of magnitude below the SN Ia rate",
    # Hollands et al. 2025 for J1637+3631
    "hollands25_vej": (1870.0, 300.0, 360.0),
    "hollands25_tflight_myr": (4.5, 0.5, 0.4),
    # Type Iax supernovae per Type Ia in a given volume (Foley et al. 2013)
    "foley13_iax_per_ia": (0.31, 0.13, 0.17),
    # Local escape speed (Deason et al. 2019)
    "vesc_sun": 528.0,
}
