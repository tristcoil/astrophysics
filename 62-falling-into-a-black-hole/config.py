"""Project configuration: paths, constants, the black holes we fall into, and rendering choices.

Every number the project *assumes* lives here with its source, so any figure or
printed result traces back to one line. Masses come from earlier projects in
this atlas when their results are on disk (H1, I4, I5), with the values they
published in their READMEs as fallbacks.

Units: the physics modules work in geometric units G = c = 1 with lengths in
units of M (= GM/c^2) unless a name says otherwise. This module converts.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
GAIA_DIR = RAW_DIR / "gaia"
CHUNK_DIR = GAIA_DIR / "chunks"
REF_DIR = RAW_DIR / "reference"
CACHE_DIR = DATA_DIR / "cache"
FRAME_DIR = CACHE_DIR / "frames"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, GAIA_DIR, CHUNK_DIR, REF_DIR, CACHE_DIR, FRAME_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: Earlier projects whose measured masses we reuse.
H1_RESULTS = ATLAS_DIR / "17-galactic-centre" / "data" / "results.json"
I4_RESULTS = ATLAS_DIR / "35-x-ray-binaries" / "data" / "results.json"
I5_RESULTS = ATLAS_DIR / "36-m87-shadow" / "data" / "results.json"

# ---------------------------------------------------------------------------
# Physical constants (SI; CODATA 2018, IAU 2015 Resolution B3)
# ---------------------------------------------------------------------------

C_LIGHT = 299_792_458.0            # m/s
GM_SUN = 1.32712440041e20          # m^3/s^2 (heliocentric gravitational constant, DE440)
G_NEWTON = 6.67430e-11             # m^3 kg^-1 s^-2 (only for kg conversions)
M_SUN_KG = GM_SUN / G_NEWTON
G0 = 9.80665                       # m/s^2, standard gravity
AU = 1.495978707e11                # m
PC = 3.0856775814913673e16         # m
YEAR = 365.25 * 86400.0            # s (Julian year)
R_SUN = 6.957e8                    # m (IAU nominal)

#: GM/c^3 for one solar mass: the time unit of a black hole, 4.925 microseconds.
T_SUN = GM_SUN / C_LIGHT**3
#: GM/c^2 for one solar mass: 1.477 km.
L_SUN = GM_SUN / C_LIGHT**2

#: Critical impact parameter of a Schwarzschild black hole (photon sphere at r = 3M).
B_CRIT = math.sqrt(27.0)

# ---------------------------------------------------------------------------
# The people we drop in: a 2 m tall body that survives ~10 g of stretching
# ---------------------------------------------------------------------------

BODY_LENGTH_M = 2.0
#: Rough survivable sustained acceleration difference head-to-foot. Fighter pilots
#: tolerate ~9 g briefly; we use 10 g as the line, as in the atlas spec.
SURVIVAL_G = 10.0

# ---------------------------------------------------------------------------
# The black holes (inputs). Each entry: mass in solar masses with 68% bounds,
# dimensionless spin where measured, where it is, and where the number comes from.
# ---------------------------------------------------------------------------

#: Fallbacks if the earlier projects' results files are missing on this machine.
_H1_FALLBACK = (3.9565, 0.0423)                 # 1e6 Msun, H1 README (joint Keck+VLT GR fit)
_I5_FALLBACK = (6.409, 5.683, 7.307)            # 1e9 Msun median, lo68, hi68 (I5 README)
_I4_FALLBACK = (21.186, 2.240, 2.259)           # Msun median, minus, plus (I4 README, Cyg X-1)

#: GW250114's remnant, as the LVK paper's main text gives it: "The remnant black hole has
#: a mass of Mf = 62.7 (+1.0 -1.1) Msun and a spin of chi_f = 0.68 (+0.01 -0.01)"
#: (source frame, 90% credible intervals; full inspiral-merger-ringdown analysis).
#: GWOSC's GWTC-5.0 preferred run lists 62.9 (-1.1 +1.0). I7 will measure it from the ringdown.
GW250114 = {
    "mass": (62.7, 1.1, 1.0),      # source-frame remnant mass, Msun (median, minus, plus; 90%)
    "spin": (0.68, 0.01, 0.01),
    "ref": "LVK 2025, PRL 135, 111403 (doi:10.1103/kw5g-d732; arXiv:2509.08054)",
    "gwosc_final_mass": (62.9, 1.1, 1.0),
}

#: Cyg X-1's spin from continuum fitting with M = 21.2, i = 27.51 deg, D = 2.22 kpc
#: ("a_* > 0.9985 (3 sigma)"; contested by some reflection/polarimetry work).
CYGX1_SPIN = {"lower_limit_3sigma": 0.9985,
              "ref": "Zhao et al. 2021, ApJ 908, 117 (doi:10.3847/1538-4357/abbcd6)"}

#: Gaia BH1: the nearest known black hole (480 pc). Joint astrometric + RV mass.
GAIA_BH1 = {"mass": (9.62, 0.18, 0.18), "distance_pc": 480.0,
            "ref": "El-Badry et al. 2023, MNRAS 518, 1057 (I4 reproduced 13.2 from Gaia DR3 alone)"}

#: Sgr A*: GRAVITY's value at their distance, for comparison with H1's own fit.
SGRA_GRAVITY = {"mass_1e6": (4.297, 0.012), "r0_kpc": 8.277, "ref": "GRAVITY Collaboration 2022, A&A 657, L12"}


def _load(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def black_holes() -> dict:
    """The black holes of this project, with masses from the atlas's own measurements.

    Returns a dict keyed by a short id; each value has ``name``, ``mass`` (median,
    minus, plus in solar masses), ``spin`` (or None), ``source`` and ``note``.
    """
    out = {}

    h1 = _load(H1_RESULTS)
    if h1 is not None:
        p = h1["fits"]["joint_gr"]["parameters"]["gm6"]
        m, s = p["median"], p["sd"]
        src = "H1 (17-galactic-centre): S2's orbit, joint Keck + VLT GR fit"
    else:
        m, s = _H1_FALLBACK
        src = "H1 README (results file not found)"
    out["sgra"] = {"name": "Sgr A*", "long_name": "Sagittarius A*, the Milky Way's central black hole",
                   "mass": (m * 1e6, s * 1e6, s * 1e6), "spin": None, "source": src,
                   "note": "spin unknown (EHT models allow 0 to ~0.9)"}

    i5 = _load(I5_RESULTS)
    if i5 is not None:
        q = i5["m87_mass"]["grmhd_calibrated"]["mass_1e9"]
        med, lo, hi = q["median"], q["lo68"], q["hi68"]
        src = "I5 (36-m87-shadow): EHT ring, GRMHD-calibrated"
    else:
        med, lo, hi = _I5_FALLBACK
        src = "I5 README (results file not found)"
    out["m87"] = {"name": "M87*", "long_name": "the black hole in galaxy M87 (EHT image)",
                  "mass": (med * 1e9, (med - lo) * 1e9, (hi - med) * 1e9), "spin": None, "source": src,
                  "note": "spin not measured directly"}

    i4 = _load(I4_RESULTS)
    if i4 is not None:
        q = i4["weighing"]["cygx1_check"]["m_bh"]
        med, mi, pl = q["median"], q["minus"], q["plus"]
        src = "I4 (35-x-ray-binaries): Cyg X-1's orbit (Miller-Jones et al. 2021 inputs)"
    else:
        med, mi, pl = _I4_FALLBACK
        src = "I4 README (results file not found)"
    out["cygx1"] = {"name": "Cyg X-1", "long_name": "Cygnus X-1, the first black hole ever identified",
                    "mass": (med, mi, pl), "spin": CYGX1_SPIN["lower_limit_3sigma"], "source": src,
                    "note": f"spin > {CYGX1_SPIN['lower_limit_3sigma']} ({CYGX1_SPIN['ref']}; contested)"}

    out["gw250114"] = {"name": "GW250114 remnant", "long_name": "the black hole born in the merger GW250114",
                       "mass": GW250114["mass"], "spin": GW250114["spin"][0], "source": GW250114["ref"],
                       "note": "90% credible bounds; I7 (not built yet) will measure it from the ringdown"}

    out["gaiabh1"] = {"name": "Gaia BH1", "long_name": "Gaia BH1, the nearest known black hole (1,560 light-years)",
                      "mass": GAIA_BH1["mass"], "spin": None, "source": GAIA_BH1["ref"],
                      "note": "dormant: no gas, found from its companion star's wobble"}
    return out


#: Order in tables and figures (biggest first).
BH_ORDER = ("m87", "sgra", "gw250114", "cygx1", "gaiabh1")

# ---------------------------------------------------------------------------
# Unit helpers
# ---------------------------------------------------------------------------


def gm_over_c2_m(mass_msun: float) -> float:
    """Gravitational radius GM/c^2 in metres."""
    return mass_msun * L_SUN


def gm_over_c3_s(mass_msun: float) -> float:
    """Light-crossing time of the gravitational radius, GM/c^3, in seconds."""
    return mass_msun * T_SUN


# ---------------------------------------------------------------------------
# The sky (Gaia DR3 + Hipparcos)
# ---------------------------------------------------------------------------

GAIA_SYNC_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
#: Stars brighter than this are drawn as points (and lensed one by one);
#: all fainter Gaia stars are summed into a HEALPix flux map.
G_POINT_LIMIT = 11.0
#: HEALPix level (nested, ICRS) of the summed-flux map: level 10 = 3.4 arcmin pixels.
FLUX_MAP_LEVEL = 10
#: Base HEALPix level for chunking queries (level 1 = 48 chunks); split on timeout.
QUERY_BASE_LEVEL = 1
QUERY_MAX_LEVEL = 5
#: Hipparcos stars brighter than this V replace their Gaia entries (Gaia's photometry
#: of the very brightest stars is saturated, and a few are missing from DR3).
HIP_V_LIMIT = 5.0

#: Gaia DR3 photometric zero points (Vega system; Riello et al. 2021, Table 3).
ZP_G, ZP_BP, ZP_RP = 25.6874, 25.3385, 24.7479
#: G - V as a function of BP - RP (Gaia DR3 documentation, Sect. 5.5.1, Table 5.9;
#: valid for -0.5 < BP-RP < 2.75): G - V = c0 + c1 x + c2 x^2 + c3 x^3.
G_MINUS_V = (-0.02704, 0.01424, -0.2156, 0.01426)

#: Sgr A* (Reid & Brunthaler 2004, ApJ 616, 872): RA 17h45m40.0409s, Dec -29d00m28.118s (ICRS).
SGRA_RA_DEG = 266.416837
SGRA_DEC_DEG = -29.007810

# ---------------------------------------------------------------------------
# Rendering choices (all of them are display choices, not physics)
# ---------------------------------------------------------------------------

#: Equirectangular sky texture size (Galactic coordinates), base level of the mipmap.
TEXTURE_WIDTH = 8192
#: Exposure: the V-band surface brightness (mag/arcsec^2) shown at display level 1
#: before tone mapping. ~19 makes the Milky Way's star clouds clearly visible, like
#: a long-exposure wide-field photograph rather than the naked eye.
EXPOSURE_MAG_ARCSEC2 = 20.5
MOVIE_FPS = 30
RANDOM_SEED = 20261019
