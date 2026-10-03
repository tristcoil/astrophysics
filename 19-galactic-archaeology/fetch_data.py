"""Download and cache everything the analysis needs.

    python fetch_data.py            # everything (about 0.75 GB, 30-60 minutes)
    python fetch_data.py --only gd1 # one piece: gd1, m92, ring, halo, gspspec, dust, literature

All Gaia data come from Gaia DR3 at ESA's archive, as many small synchronous
queries (see gaia_archive.py for why). Each sample is designed so the cut that
matters is made on the server and every other decision is made locally, where
it can be inspected and changed without re-downloading:

1. **The GD-1 field** (~1.7 million stars): every star brighter than G = 20.5
   in a 95 x 14 degree strip along the stream, minus obvious nearby stars
   (parallax > 1 mas) and red dwarfs (BP-RP > 1.6). Plus a server-side count of
   *all* stars in the strip, binned on the sky -- the "before" picture.
2. **M92** (a globular cluster as old and metal-poor as GD-1's parent), whose
   colour-magnitude diagram becomes the filter that picks GD-1's stars.
3. **The ring** (~2.9 million stars): every star with a Gaia radial velocity
   and a parallax better than 20%, whose Galactocentric radius is within
   0.1 kpc of the Sun's -- Antoja et al.'s (2018) phase-spiral sample, in DR3.
4. **The fast stars** (~1 million): 6D stars moving > 150 km/s relative to the
   local disc, plus a random 2% of all 6D stars for context, with Gaia's
   spectroscopic metallicities and Bailer-Jones et al. (2021) distances joined
   on the server.
5. **The metal-poor stars**: 6D stars whose Gaia RVS spectra say [M/H] < -0.7,
   with *no* cut on velocity -- the route to the halo that avoids circularity.
6. **Dust**: the Schlegel, Finkbeiner & Davis (1998) reddening map.
7. **Literature**: spectroscopically confirmed GD-1 members (Huang et al.
   2019), used only to grade the membership.

The server-side counts that each sample is checked against are recorded in a
JSON file next to it, with the exact ADQL and the retrieval date.
"""

from __future__ import annotations

import argparse
import io
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

import gaia_archive as ga
from config import (
    AUX_DIR, CHUNK_DIR, GAIA_DIR, GD1_BP_RP_MAX, GD1_COUNT_MAP_BIN_DEG, GD1_G_MAX,
    GD1_PARALLAX_MAX_MAS, GD1_PHI1_RANGE, GD1_PHI2_RANGE, GD1_QUERY_SEGMENT_DEG,
    GD1_ROTATION, GSPSPEC_FETCH_MH_MAX, HALO_FETCH_SPEED_KMS, K_KMS, M92,
    N_DR3_SOURCES, PARALLAX_OVER_ERROR_6D, R0_KPC, RANDOM_FRACTION, RING_HALF_WIDTH_KPC,
    SFD_URL, V_PECULIAR_KMS,
)

# --------------------------------------------------------------------------
# ADQL building blocks
# --------------------------------------------------------------------------

#: ICRS -> Galactic rotation matrix (astropy's, to 16 digits).
ICRS_TO_GAL = (
    (-0.0548756577125916, -0.8734370519556159, -0.4838350736167155),
    (0.4941094371927268, -0.4448297212232952, 0.7469821839866676),
    (-0.8676661375596576, -0.1980763372730005, 0.4559838136873016),
)


def _f(x: float) -> str:
    return repr(float(x))


def adql_unit_vector(p: str = "") -> tuple[str, str, str]:
    """ICRS unit vector of a star, as three ADQL expressions."""
    ra, de = f"RADIANS({p}ra)", f"RADIANS({p}dec)"
    return (f"COS({de})*COS({ra})", f"COS({de})*SIN({ra})", f"SIN({de})")


def adql_gd1(p: str = "") -> tuple[str, str]:
    """GD-1 stream coordinates (phi1, phi2) in degrees, as ADQL."""
    u = adql_unit_vector(p)
    rows = [f"({' + '.join(f'{_f(R[i])}*{u[i]}' for i in range(3))})" for R in GD1_ROTATION]
    return f"DEGREES(ATAN2({rows[1]}, {rows[0]}))", f"DEGREES(ASIN({rows[2]}))"


def adql_uvw(p: str = "") -> tuple[str, str, str]:
    """Heliocentric Galactic velocity (U towards the Galactic Centre, V along
    rotation, W towards the North Galactic Pole), km/s, as ADQL. Distance is
    1/parallax. Verified against astropy to 1e-6 km/s (test_physics.py)."""
    ra, de = f"RADIANS({p}ra)", f"RADIANS({p}dec)"
    va, vd = f"({_f(K_KMS)}*{p}pmra/{p}parallax)", f"({_f(K_KMS)}*{p}pmdec/{p}parallax)"
    rv = f"{p}radial_velocity"
    vx = f"({rv}*COS({de})*COS({ra}) - {va}*SIN({ra}) - {vd}*SIN({de})*COS({ra}))"
    vy = f"({rv}*COS({de})*SIN({ra}) + {va}*COS({ra}) - {vd}*SIN({de})*SIN({ra}))"
    vz = f"({rv}*SIN({de}) + {vd}*COS({de}))"
    return tuple(
        f"({_f(T[0])}*{vx} + {_f(T[1])}*{vy} + {_f(T[2])}*{vz})" for T in ICRS_TO_GAL
    )


def adql_speed_vs_lsr(p: str = "") -> str:
    """|v - v_LSR|: speed relative to a circular orbit at the Sun, km/s."""
    u, v, w = adql_uvw(p)
    pu, pv, pw = V_PECULIAR_KMS
    return (f"SQRT(POWER({u} + {_f(pu)}, 2) + POWER({v} + {_f(pv)}, 2) "
            f"+ POWER({w} + {_f(pw)}, 2))")


def adql_galactocentric_radius(p: str = "") -> str:
    """Galactocentric cylindrical radius in kpc, with distance = 1/parallax and
    the Sun in the midplane (its 21 pc height changes R by < 3 pc here)."""
    x = f"COS(RADIANS({p}b))*COS(RADIANS({p}l))/{p}parallax"
    y = f"COS(RADIANS({p}b))*SIN(RADIANS({p}l))/{p}parallax"
    return f"SQRT(POWER({_f(R0_KPC)} - {x}, 2) + POWER({y}, 2))"


def _gd1_to_icrs(phi1: float, phi2: float) -> tuple[float, float]:
    R = np.array(GD1_ROTATION)
    p1, p2 = np.radians(phi1), np.radians(phi2)
    v = R.T @ np.array([np.cos(p1) * np.cos(p2), np.sin(p1) * np.cos(p2), np.sin(p2)])
    return float(np.degrees(np.arctan2(v[1], v[0])) % 360), float(np.degrees(np.arcsin(v[2])))


def _gd1_segment_circle(p1a: float, p1b: float, margin: float = 0.3) -> str:
    """Smallest sky circle enclosing a phi1 segment of the field (+ margin)."""
    p2a, p2b = GD1_PHI2_RANGE
    ra, dec = _gd1_to_icrs(0.5 * (p1a + p1b), 0.5 * (p2a + p2b))
    radius = np.hypot(0.5 * (p1b - p1a), 0.5 * (p2b - p2a)) + margin
    return f"CIRCLE('ICRS', {ra:.5f}, {dec:.5f}, {radius:.3f})"


def _gd1_segments() -> list[tuple[float, float]]:
    lo, hi = GD1_PHI1_RANGE
    edges = np.arange(lo, hi + 1e-9, GD1_QUERY_SEGMENT_DEG)
    return list(zip(edges[:-1], edges[1:]))


# --------------------------------------------------------------------------
# Column lists
# --------------------------------------------------------------------------

ASTROMETRY = [
    "source_id", "ra", "dec", "parallax", "parallax_error", "pmra", "pmra_error",
    "pmdec", "pmdec_error", "ruwe", "astrometric_params_solved",
    "nu_eff_used_in_astrometry", "pseudocolour", "phot_g_mean_mag", "bp_rp",
]
RV = ["radial_velocity", "radial_velocity_error"]
GD1_COLUMNS = ASTROMETRY + [
    "pmra_pmdec_corr", "phot_bp_mean_mag", "phot_rp_mean_mag",
    "phot_bp_rp_excess_factor", "visibility_periods_used",
] + RV
M92_COLUMNS = GD1_COLUMNS
RING_COLUMNS = ASTROMETRY + RV
GSPSPEC = [
    "mh_gspspec", "mh_gspspec_lower", "mh_gspspec_upper", "alphafe_gspspec",
    "alphafe_gspspec_lower", "alphafe_gspspec_upper", "fem_gspspec", "teff_gspspec",
    "logg_gspspec", "flags_gspspec",
]
BAILER_JONES = [
    "r_med_geo", "r_lo_geo", "r_hi_geo", "r_med_photogeo", "r_lo_photogeo", "r_hi_photogeo",
]
HALO_COLUMNS = (
    [f"g.{c}" for c in ASTROMETRY + RV + ["phot_bp_rp_excess_factor", "random_index"]]
    + [f"ap.{c}" for c in GSPSPEC] + [f"d.{c}" for c in BAILER_JONES]
)

# --------------------------------------------------------------------------
# 1. The GD-1 field
# --------------------------------------------------------------------------

GD1_PATH = GAIA_DIR / "gd1_field.fits"
GD1_COUNTS_PATH = GAIA_DIR / "gd1_countmap.fits"


def gd1_query(p1a: float, p1b: float) -> str:
    phi1, phi2 = adql_gd1()
    return (
        f"SELECT {', '.join(GD1_COLUMNS)} FROM gaiadr3.gaia_source\n"
        f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), {_gd1_segment_circle(p1a, p1b)})\n"
        f"  AND phot_g_mean_mag < {GD1_G_MAX} AND parallax < {GD1_PARALLAX_MAX_MAS}"
        f" AND bp_rp < {GD1_BP_RP_MAX}\n"
        f"  AND {phi1} >= {p1a} AND {phi1} < {p1b}\n"
        f"  AND {phi2} BETWEEN {GD1_PHI2_RANGE[0]} AND {GD1_PHI2_RANGE[1]}"
    )


def gd1_count_query(p1a: float, p1b: float) -> str:
    """All stars (no cuts but G) binned on the sky in stream coordinates."""
    phi1, phi2 = adql_gd1()
    b = GD1_COUNT_MAP_BIN_DEG
    return (
        f"SELECT FLOOR({phi1}/{b}) AS i1, FLOOR({phi2}/{b}) AS i2, COUNT(*) AS n\n"
        f"FROM gaiadr3.gaia_source\n"
        f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), {_gd1_segment_circle(p1a, p1b)})\n"
        f"  AND phot_g_mean_mag < {GD1_G_MAX}\n"
        f"  AND {phi1} >= {p1a} AND {phi1} < {p1b}\n"
        f"  AND {phi2} BETWEEN {GD1_PHI2_RANGE[0]} AND {GD1_PHI2_RANGE[1]}\n"
        f"GROUP BY i1, i2"
    )


def _fetch_segments(name: str, make_query) -> tuple[Table, list[dict]]:
    tables, log = [], []
    for p1a, p1b in _gd1_segments():
        path = CHUNK_DIR / name / f"phi1_{p1a:+06.1f}.fits"
        if path.exists():
            t = Table.read(path, format="fits")
        else:
            t0 = time.time()
            t = ga.sync_query(make_query(p1a, p1b), timeout=240)
            path.parent.mkdir(parents=True, exist_ok=True)
            t.write(path, format="fits", overwrite=True)
            print(f"    phi1 {p1a:+6.1f} .. {p1b:+6.1f}: {len(t):>9,} rows in "
                  f"{time.time() - t0:5.1f} s", flush=True)
        tables.append(t)
        log.append({"phi1": [p1a, p1b], "rows": len(t)})
    return vstack(tables, metadata_conflicts="silent"), log


def fetch_gd1(refresh: bool = False) -> None:
    if not GD1_PATH.exists() or refresh:
        print("GD-1 field: stars along the stream")
        table, log = _fetch_segments("gd1_field", gd1_query)
        if len(np.unique(table["source_id"])) != len(table):
            raise RuntimeError("GD-1 field: duplicate source_ids across segments")
        ga.write_with_provenance(table, GD1_PATH, {
            "what": "Gaia DR3 stars in the GD-1 field",
            "service": "ESA Gaia archive, synchronous TAP",
            "retrieved_utc": ga.now_utc(),
            "example_query": gd1_query(*_gd1_segments()[0]),
            "segments": log,
        })
        _drop_chunks("gd1_field")
    if not GD1_COUNTS_PATH.exists() or refresh:
        print("GD-1 field: server-side star counts (the unfiltered sky)")
        table, log = _fetch_segments("gd1_countmap", gd1_count_query)
        ga.write_with_provenance(table, GD1_COUNTS_PATH, {
            "what": f"counts of all Gaia DR3 sources with G < {GD1_G_MAX} in "
                    f"{GD1_COUNT_MAP_BIN_DEG}-degree bins of (phi1, phi2); bin index "
                    f"i = floor(phi / {GD1_COUNT_MAP_BIN_DEG})",
            "retrieved_utc": ga.now_utc(),
            "example_query": gd1_count_query(*_gd1_segments()[0]),
            "segments": log,
            "total_stars": int(np.sum(table["n"])),
        })
        _drop_chunks("gd1_countmap")


# --------------------------------------------------------------------------
# 2. M92
# --------------------------------------------------------------------------

M92_PATH = GAIA_DIR / "m92.fits"


def fetch_m92(refresh: bool = False) -> None:
    if M92_PATH.exists() and not refresh:
        return
    print("M92: the template cluster")
    q = (f"SELECT {', '.join(M92_COLUMNS)} FROM gaiadr3.gaia_source\n"
         f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), "
         f"CIRCLE('ICRS', {M92.ra_deg}, {M92.dec_deg}, {M92.radius_deg}))")
    table = ga.sync_query(q)
    ga.write_with_provenance(table, M92_PATH, {
        "what": f"all Gaia DR3 sources within {M92.radius_deg} deg of {M92.name}",
        "retrieved_utc": ga.now_utc(), "query": q,
    })
    print(f"    {len(table):,} stars")


# --------------------------------------------------------------------------
# 3-5. The all-sky 6D samples
# --------------------------------------------------------------------------

RING_PATH = GAIA_DIR / "ring.fits"
HALO_PATH = GAIA_DIR / "fast_and_random.fits"
GSPSPEC_PATH = GAIA_DIR / "metal_poor.fits"

#: Fetched slightly wider than used: the server's R neglects the Sun's 21 pc
#: height, which moves R by up to ~3 pc for stars 1 kpc from the plane.
_RING_MARGIN_KPC = 0.005

SIX_D = (f"radial_velocity IS NOT NULL AND parallax_over_error > {PARALLAX_OVER_ERROR_6D}")


def ring_template(select: str) -> str:
    r = adql_galactocentric_radius()
    lo, hi = R0_KPC - RING_HALF_WIDTH_KPC - _RING_MARGIN_KPC, R0_KPC + RING_HALF_WIDTH_KPC + _RING_MARGIN_KPC
    return (f"SELECT {select} FROM gaiadr3.gaia_source\n"
            f"WHERE source_id >= {{lo}} AND source_id < {{hi}} AND {SIX_D}\n"
            f"  AND {r} BETWEEN {lo:.4f} AND {hi:.4f}")


def halo_template(select: str) -> str:
    random_cut = int(round(RANDOM_FRACTION * N_DR3_SOURCES))
    joins = ("" if select.startswith("COUNT") else
             "LEFT OUTER JOIN gaiadr3.astrophysical_parameters AS ap ON ap.source_id = g.source_id\n"
             "LEFT OUTER JOIN external.gaiaedr3_distance AS d ON d.source_id = g.source_id\n")
    return (f"SELECT {select} FROM gaiadr3.gaia_source AS g\n{joins}"
            f"WHERE g.source_id >= {{lo}} AND g.source_id < {{hi}}\n"
            f"  AND g.radial_velocity IS NOT NULL AND g.parallax_over_error > {PARALLAX_OVER_ERROR_6D}\n"
            f"  AND ({adql_speed_vs_lsr('g.')} > {HALO_FETCH_SPEED_KMS}"
            f" OR g.random_index < {random_cut})")


def gspspec_template(select: str) -> str:
    bj = ("" if select.startswith("COUNT") else
          "LEFT OUTER JOIN external.gaiaedr3_distance AS d ON d.source_id = g.source_id\n")
    return (f"SELECT {select} FROM gaiadr3.gaia_source AS g\n"
            f"JOIN gaiadr3.astrophysical_parameters AS ap ON ap.source_id = g.source_id\n{bj}"
            f"WHERE g.source_id >= {{lo}} AND g.source_id < {{hi}}"
            f" AND ap.source_id >= {{lo}} AND ap.source_id < {{hi}}\n"
            f"  AND g.radial_velocity IS NOT NULL AND g.parallax_over_error > {PARALLAX_OVER_ERROR_6D}\n"
            f"  AND ap.mh_gspspec < {GSPSPEC_FETCH_MH_MAX}")


SAMPLES = {
    "ring": (RING_PATH, ring_template, ", ".join(RING_COLUMNS),
             "6D stars within 0.1 kpc of the Sun's Galactocentric radius"),
    "halo": (HALO_PATH, halo_template, ", ".join(HALO_COLUMNS),
             f"6D stars faster than {HALO_FETCH_SPEED_KMS:.0f} km/s relative to the LSR, "
             f"plus a random {RANDOM_FRACTION:.0%} of all 6D stars"),
    "gspspec": (GSPSPEC_PATH, gspspec_template,
                ", ".join([f"g.{c}" for c in ASTROMETRY + RV + ["phot_bp_rp_excess_factor",
                                                                   "random_index"]]
                          + [f"ap.{c}" for c in GSPSPEC] + [f"d.{c}" for c in BAILER_JONES]),
                f"6D stars with GSP-Spec [M/H] < {GSPSPEC_FETCH_MH_MAX} (no kinematic cut)"),
}


def fetch_sample(key: str, refresh: bool = False, verify: bool = True) -> None:
    path, template, select, what = SAMPLES[key]
    if path.exists() and not refresh:
        return
    print(f"{key}: {what}")
    t0 = time.time()
    table, log = ga.fetch_all_sky(key, template(select))
    if len(np.unique(table["source_id"])) != len(table):
        raise RuntimeError(f"{key}: duplicate source_ids")
    provenance = {
        "what": what, "service": "ESA Gaia archive, synchronous TAP",
        "retrieved_utc": ga.now_utc(), "query_template": template(select),
        "download_minutes": round((time.time() - t0) / 60, 1),
        "chunks": [vars(c) for c in log],
    }
    if verify:
        print(f"  verifying against a server-side count ...", flush=True)
        n_server = ga.count_all_sky(template("COUNT(*)"))
        provenance["server_count"] = n_server
        if n_server != len(table):
            raise RuntimeError(f"{key}: {len(table):,} rows downloaded but the archive "
                               f"counts {n_server:,}")
        print(f"  {len(table):,} rows = server count", flush=True)
    ga.write_with_provenance(table, path, provenance)
    _drop_chunks(key)


def _drop_chunks(name: str) -> None:
    """The per-chunk cache exists only so an interrupted download can resume;
    once the merged table is verified and written it is redundant."""
    import shutil

    shutil.rmtree(CHUNK_DIR / name, ignore_errors=True)


# --------------------------------------------------------------------------
# 6. Dust
# --------------------------------------------------------------------------

SFD_PATH = AUX_DIR / "lambda_sfd_ebv.fits"


def _download(url: str, path: Path, *, timeout: float = 600) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (astro-atlas H3)"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    tmp = path.with_suffix(".part")
    tmp.write_bytes(data)
    tmp.rename(path)


def fetch_dust(refresh: bool = False) -> None:
    if SFD_PATH.exists() and not refresh:
        return
    print("dust: SFD (1998) E(B-V), HEALPix, from NASA LAMBDA")
    _download(SFD_URL, SFD_PATH)
    SFD_PATH.with_suffix(".json").write_text(json.dumps(
        {"url": SFD_URL, "retrieved_utc": ga.now_utc(),
         "bytes": SFD_PATH.stat().st_size}, indent=2) + "\n")


# --------------------------------------------------------------------------
# 7. Literature
# --------------------------------------------------------------------------

HUANG_PATH = AUX_DIR / "huang2019_gd1_members.tsv"
VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"


def fetch_literature(refresh: bool = False) -> None:
    if HUANG_PATH.exists() and not refresh:
        return
    print("literature: Huang et al. (2019) spectroscopic GD-1 members (VizieR J/ApJ/877/13)")
    url = VIZIER_ASU + "?" + urllib.parse.urlencode(
        {"-source": "J/ApJ/877/13", "-out.all": "", "-out.max": "unlimited"})
    _download(url, HUANG_PATH)


# --------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", action="append",
                        choices=["gd1", "m92", "ring", "halo", "gspspec", "dust", "literature"])
    parser.add_argument("--refresh", action="store_true", help="re-download even if cached")
    parser.add_argument("--no-verify", action="store_true",
                        help="skip the server-side row counts")
    args = parser.parse_args(argv)
    todo = args.only or ["dust", "literature", "m92", "gd1", "gspspec", "halo", "ring"]
    for key in todo:
        if key == "gd1":
            fetch_gd1(args.refresh)
        elif key == "m92":
            fetch_m92(args.refresh)
        elif key == "dust":
            fetch_dust(args.refresh)
        elif key == "literature":
            fetch_literature(args.refresh)
        else:
            fetch_sample(key, args.refresh, verify=not args.no_verify)
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
