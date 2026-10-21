"""Download and cache everything the analysis needs.

Seven datasets, all free and public:

1. **Gaia DR3** (ESA): every source with parallax > 10 mas measured to better
   than 10%, joined with the Rybizki et al. (2022) astrometric-fidelity
   classifier and flagged against the Gaia Catalogue of Nearby Stars.
2. **The Gaia Catalogue of Nearby Stars** (GCNS; Gaia Collaboration, Smart et al.
   2021), the Gaia team's own curated 100 pc census, used only to check our cuts.
3. **Hipparcos** (ESA, 1989-93; van Leeuwen 2007 reduction). Gaia cannot measure
   the very brightest stars -- its detectors saturate -- so Sirius, Vega, Alpha
   Centauri and friends have to come from its predecessor.
4. **The 10 parsec census** (Reyle et al. 2021, updated 2023): the best current
   inventory of every star and brown dwarf within 10 pc, from all sources. It
   tells us what Gaia misses.
5. **Mamajek's dwarf sequence**: the standard table linking colour, temperature,
   spectral type and luminosity for main-sequence stars.
6. **SIMBAD** (CDS): names for nearby stars, and positions for famous ones.
7. **NASA Exoplanet Archive**: which of these stars have known planets.

Everything is cached under ``data/``, with the retrieval date and the exact query
recorded in ``data/provenance.json``. Nothing is re-downloaded that is already
present, and the analysis runs offline once the cache is warm.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import warnings
from datetime import datetime, timezone

import numpy as np
from astropy.table import Table, vstack

import config as C

warnings.filterwarnings("ignore", category=UserWarning)

GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
DWARF_SEQUENCE_URL = "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt"
CENSUS_10PC_VIZIER = "J/A+A/650/A201/tablea1"

#: Gaia source_ids encode sky position: the top bits are a HEALPix index, so
#: ``source_id // 2**57`` is the level-1 HEALPix pixel (0..47). Splitting the
#: query on those boundaries cuts the sky into 48 equal-area patches of a few
#: thousand stars each -- small enough that each is a quick synchronous query.
N_SKY_CHUNKS = 48
CHUNK_STEP = 2**57

GAIA_COLUMNS = [
    "g.source_id", "g.ra", "g.dec", "g.l", "g.b",
    "g.parallax", "g.parallax_error",
    "g.pmra", "g.pmra_error", "g.pmdec", "g.pmdec_error",
    "g.radial_velocity", "g.radial_velocity_error",
    "g.phot_g_mean_mag", "g.phot_bp_mean_mag", "g.phot_rp_mean_mag", "g.bp_rp",
    "g.phot_bp_rp_excess_factor",
    "g.ruwe", "g.astrometric_params_solved", "g.non_single_star",
    "s.fidelity_v2",
    "c.gcns_prob", "c.wd_prob",
]

#: Where-clause shared by the data query and the count that verifies it.
GAIA_WHERE = (f"g.parallax > {C.MIN_PARALLAX_MAS:g} "
              f"AND g.parallax_over_error > {C.MIN_PARALLAX_OVER_ERROR:g}")

GCNS_COLUMNS = [
    "source_id", "ra", "dec", "parallax", "parallax_error", "phot_g_mean_mag",
    "phot_bp_mean_mag", "phot_rp_mean_mag", "dist_16", "dist_50", "dist_84",
    "gcns_prob", "wd_prob",
]

HIPPARCOS_QUERY = f"""
SELECT n.hip, n.ra, n.dec, n.plx, n.e_plx, n.pm_ra, n.pm_de, n.e_pm_ra, n.e_pm_de,
       n.b_v, n.e_b_v, n.v_i, n.hp_mag, n.sn, n.nc,
       h.vmag, h.sptype, h.multflag, h.ccdm,
       x.source_id AS gaia_source_id, x.angular_distance AS xm_sep_arcsec,
       g.parallax AS gaia_parallax, g.parallax_over_error AS gaia_poe,
       g.phot_g_mean_mag AS gaia_g
FROM public.hipparcos_newreduction AS n
LEFT OUTER JOIN public.hipparcos AS h ON h.hip = n.hip
LEFT OUTER JOIN gaiadr3.hipparcos2_best_neighbour AS x ON x.original_ext_source_id = n.hip
LEFT OUTER JOIN gaiadr3.gaia_source AS g ON g.source_id = x.source_id
WHERE n.plx > {C.MIN_PARALLAX_MAS:g}
"""

#: Every star in the sky bright enough to see without a telescope, whatever its
#: distance: the answer to "how many of the stars I can see are on this map?"
NAKED_EYE_QUERY = """
SELECT h.hip, h.vmag, h.b_v, h.sptype, n.plx, n.e_plx
FROM public.hipparcos AS h
LEFT OUTER JOIN public.hipparcos_newreduction AS n ON n.hip = h.hip
WHERE h.vmag < 6.5
"""

#: Columns pulled from the NASA Exoplanet Archive's composite table.
EXOPLANET_COLUMNS = [
    "pl_name", "hostname", "ra", "dec", "sy_dist", "sy_pnum", "sy_snum",
    "discoverymethod", "disc_year", "pl_rade", "pl_bmasse", "pl_orbper",
    "st_spectype", "gaia_dr3_id", "hip_name",
]


# --------------------------------------------------------------------------
# Plumbing
# --------------------------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _record_provenance(key: str, **info) -> None:
    """Append one dataset's retrieval record to data/provenance.json."""
    prov = {}
    if C.PROVENANCE_PATH.exists():
        prov = json.loads(C.PROVENANCE_PATH.read_text())
    prov[key] = {"retrieved_utc": _now(), **info}
    C.PROVENANCE_PATH.write_text(json.dumps(prov, indent=2))


def load_provenance() -> dict:
    if not C.PROVENANCE_PATH.exists():
        return {}
    return json.loads(C.PROVENANCE_PATH.read_text())


def _lowercase(table: Table) -> Table:
    """Lower-case every column name.

    The archive is not consistent about this: the same query, sent to the same
    endpoint, came back with its first column named ``SOURCE_ID`` from some
    requests and ``source_id`` from others (presumably different servers behind
    the load balancer). Stacked naively, that splits one column into two
    half-empty ones -- every row present, so a row count does not notice.
    """
    for name in table.colnames:
        if name != name.lower():
            table.rename_column(name, name.lower())
    return table


def _fix_strings(table: Table) -> Table:
    """Turn variable-length string columns into fixed-width ones FITS can store."""
    for name in table.colnames:
        if table[name].dtype == object:
            table[name] = [str(v) if v is not None and v is not np.ma.masked else ""
                           for v in table[name]]
    return table


def gaia_sync(adql: str, *, attempts: int = 4, timeout: int = 900) -> Table:
    """Run a synchronous ADQL query against the ESA Gaia archive.

    Synchronous rather than asynchronous because, while this was written, the
    archive's asynchronous job service was failing (it is being reworked for
    Gaia DR4); plain synchronous requests were reliable. Two protections make
    that safe:

    * The archive marks a result it had to cut short with QUERY_STATUS =
      OVERFLOW. We refuse any such result rather than silently losing rows.
    * Queries are kept small (sky chunks), and the chunk total is checked
      against a server-side COUNT(*) afterwards.
    """
    data = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "QUERY": adql}
    ).encode()
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(GAIA_TAP_SYNC, data=data, timeout=timeout) as r:
                raw = r.read()
            if b'value="OVERFLOW"' in raw:
                raise RuntimeError("Gaia archive truncated the result (OVERFLOW)")
            return _lowercase(Table.read(io.BytesIO(raw), format="votable"))
        except urllib.error.HTTPError as exc:
            body = exc.read().decode(errors="replace")
            found = re.findall(r'<INFO name="QUERY_STATUS" value="ERROR">(.*?)</INFO>', body, re.S)
            message = (found[0] if found else body[:300]).strip()
            if exc.code < 500:  # a malformed query will not fix itself
                raise RuntimeError(f"Gaia archive rejected the query: {message}") from None
            error = f"HTTP {exc.code}: {message[:120]}"
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            error = f"{type(exc).__name__}: {exc}"
        if attempt == attempts:
            raise RuntimeError(f"Gaia archive query failed after {attempts} attempts: {error}")
        wait = 20 * attempt
        print(f"    archive hiccup ({error}); retrying in {wait} s", flush=True)
        time.sleep(wait)
    raise AssertionError("unreachable")


def gaia_count(from_where: str) -> int:
    """COUNT(*) over a FROM ... WHERE ... clause, computed by the archive."""
    t = gaia_sync(f"SELECT COUNT(*) AS n {from_where}")
    return int(t["n"][0])


def _chunked_download(label: str, chunk_dir, make_query, *, refresh: bool) -> Table:
    """Download a query sky-chunk by sky-chunk, caching each chunk to disk."""
    tables = []
    t0 = time.time()
    for k in range(N_SKY_CHUNKS):
        path = chunk_dir / f"chunk_{k:02d}.fits"
        if path.exists() and not refresh:
            tables.append(_lowercase(Table.read(path)))
            continue
        lo, hi = k * CHUNK_STEP, (k + 1) * CHUNK_STEP
        table = _fix_strings(gaia_sync(make_query(lo, hi)))
        table.write(path, overwrite=True)
        tables.append(table)
        print(f"    {label} chunk {k + 1:2d}/{N_SKY_CHUNKS}: {len(table):6,} rows "
              f"({time.time() - t0:.0f} s elapsed)", flush=True)
    names = {tuple(t.colnames) for t in tables}
    if len(names) != 1:
        raise RuntimeError(f"{label} chunks disagree on their columns: {names}")
    return vstack(tables, metadata_conflicts="silent")


def _clear_chunks(chunk_dir) -> None:
    """Delete per-chunk files once the combined, verified table is on disk."""
    for path in chunk_dir.glob("chunk_*.fits"):
        path.unlink()


# --------------------------------------------------------------------------
# 1. Gaia DR3 within 100 pc
# --------------------------------------------------------------------------

def _gaia_query(lo: int, hi: int) -> str:
    return (
        f"SELECT {', '.join(GAIA_COLUMNS)}\n"
        "FROM gaiadr3.gaia_source AS g\n"
        "LEFT OUTER JOIN external.gaiaedr3_spurious AS s ON s.source_id = g.source_id\n"
        "LEFT OUTER JOIN external.gaiaedr3_gcns_main_1 AS c ON c.source_id = g.source_id\n"
        f"WHERE g.source_id >= {lo} AND g.source_id < {hi}\n  AND {GAIA_WHERE}"
    )


def fetch_gaia(*, refresh: bool = False) -> Table:
    """Every Gaia DR3 source with parallax > 10 mas at better than 10% precision.

    The fidelity classifier and the GCNS membership come along as joined
    columns, so the spurious-source cut is applied (and reported) locally rather
    than hidden inside the query.

    EDR3 and DR3 share source identifiers and astrometry, so tables built on
    EDR3 (the fidelity classifier, GCNS) join directly onto DR3.
    """
    if C.GAIA_PATH.exists() and not refresh:
        return Table.read(C.GAIA_PATH)
    print("  Gaia DR3: parallax > 10 mas, parallax_over_error > 10, in 48 sky chunks")
    table = _chunked_download("Gaia", C.GAIA_CHUNK_DIR, _gaia_query, refresh=refresh)

    expected = gaia_count(f"FROM gaiadr3.gaia_source AS g WHERE {GAIA_WHERE}")
    if len(table) != expected:
        raise RuntimeError(f"Gaia download incomplete: {len(table)} rows, archive "
                           f"counts {expected}. Delete data/gaia_chunks/ and retry.")
    if np.ma.getmaskarray(np.ma.asarray(table["source_id"])).any():
        raise RuntimeError("Gaia download has rows without a source_id")
    print(f"  {len(table):,} rows, matching the archive's own count")
    table.meta.clear()
    table.write(C.GAIA_PATH, overwrite=True)
    _clear_chunks(C.GAIA_CHUNK_DIR)
    _record_provenance(
        "gaia_dr3", source="ESA Gaia archive, gaiadr3.gaia_source joined with "
        "external.gaiaedr3_spurious (Rybizki+2022) and external.gaiaedr3_gcns_main_1",
        query=_gaia_query(0, CHUNK_STEP).replace(f"g.source_id < {CHUNK_STEP}", "g.source_id < <chunk end>"),
        rows=len(table), verified_against_count=expected,
    )
    return table


# --------------------------------------------------------------------------
# 2. The Gaia Catalogue of Nearby Stars
# --------------------------------------------------------------------------

def _gcns_query(lo: int, hi: int) -> str:
    return (f"SELECT {', '.join(GCNS_COLUMNS)}\nFROM external.gaiaedr3_gcns_main_1\n"
            f"WHERE source_id >= {lo} AND source_id < {hi}")


def fetch_gcns(*, refresh: bool = False) -> Table:
    """The GCNS main table: 331,312 objects the Gaia team places within ~100 pc."""
    if C.GCNS_PATH.exists() and not refresh:
        return Table.read(C.GCNS_PATH)
    print("  GCNS: the Gaia Collaboration's curated 100 pc catalogue, in 48 sky chunks")
    table = _chunked_download("GCNS", C.GCNS_CHUNK_DIR, _gcns_query, refresh=refresh)
    expected = gaia_count("FROM external.gaiaedr3_gcns_main_1")
    if len(table) != expected:
        raise RuntimeError(f"GCNS download incomplete: {len(table)} of {expected} rows")
    if np.ma.getmaskarray(np.ma.asarray(table["source_id"])).any():
        raise RuntimeError("GCNS download has rows without a source_id")
    print(f"  {len(table):,} rows, matching the archive's own count")
    table.meta.clear()
    table.write(C.GCNS_PATH, overwrite=True)
    _clear_chunks(C.GCNS_CHUNK_DIR)
    _record_provenance("gcns", source="Gaia Collaboration, Smart et al. 2021, A&A 649, A6; "
                       "ESA Gaia archive table external.gaiaedr3_gcns_main_1",
                       rows=len(table))
    return table


# --------------------------------------------------------------------------
# 3. Hipparcos, for the stars Gaia is too dazzled to measure
# --------------------------------------------------------------------------

def fetch_hipparcos(*, refresh: bool = False) -> Table:
    """Hipparcos stars with parallax > 10 mas, with their best Gaia counterpart.

    Parallaxes come from van Leeuwen's 2007 re-reduction, which improved on the
    1997 catalogue; V magnitudes and spectral types only exist in the original,
    so both are joined. The Gaia columns come from the archive's official
    Hipparcos-2 cross-match, and tell the analysis whether Gaia already has
    this star.
    """
    if C.HIPPARCOS_PATH.exists() and not refresh:
        return Table.read(C.HIPPARCOS_PATH)
    print("  Hipparcos: parallax > 10 mas, with the Gaia DR3 cross-match")
    table = _fix_strings(gaia_sync(HIPPARCOS_QUERY))
    print(f"  {len(table):,} Hipparcos stars")
    table.meta.clear()
    table.write(C.HIPPARCOS_PATH, overwrite=True)
    _record_provenance("hipparcos", source="Hipparcos new reduction (van Leeuwen 2007) + "
                       "original catalogue (ESA 1997) + gaiadr3.hipparcos2_best_neighbour, "
                       "via the ESA Gaia archive", query=HIPPARCOS_QUERY.strip(), rows=len(table))
    return table


def fetch_naked_eye(*, refresh: bool = False) -> Table:
    """Every Hipparcos star with V < 6.5 -- the whole naked-eye sky, near or far."""
    if C.NAKED_EYE_PATH.exists() and not refresh:
        return Table.read(C.NAKED_EYE_PATH)
    print("  Hipparcos: every naked-eye star (V < 6.5), at any distance")
    table = _fix_strings(gaia_sync(NAKED_EYE_QUERY))
    print(f"  {len(table):,} stars")
    table.meta.clear()
    table.write(C.NAKED_EYE_PATH, overwrite=True)
    _record_provenance("naked_eye", source="Hipparcos catalogue (ESA 1997) V magnitudes with "
                       "van Leeuwen (2007) parallaxes, via the ESA Gaia archive",
                       query=NAKED_EYE_QUERY.strip(), rows=len(table))
    return table


# --------------------------------------------------------------------------
# 4. The 10 parsec census
# --------------------------------------------------------------------------

def fetch_census_10pc(*, refresh: bool = False) -> Table:
    """Reyle et al. (2021) 10 pc catalogue, in its 2023 update, from VizieR.

    It combines Gaia with ground-based parallaxes and infrared surveys, so it
    includes the brown dwarfs and bright stars that Gaia alone misses.
    """
    if C.CENSUS_10PC_PATH.exists() and not refresh:
        return Table.read(C.CENSUS_10PC_PATH)
    from astroquery.vizier import Vizier

    print("  10 pc census: Reyle et al. 2021 (2023 update), via VizieR")
    table = Vizier(row_limit=-1, columns=["**"]).get_catalogs(CENSUS_10PC_VIZIER)[0]
    description = table.meta.get("description", "")
    table = _fix_strings(Table(table, masked=True))
    table.meta.clear()
    table.write(C.CENSUS_10PC_PATH, overwrite=True)
    print(f"  {len(table)} entries ({description})")
    _record_provenance("census_10pc", source=f"Reyle et al. 2021, A&A 650, A201; VizieR "
                       f"{CENSUS_10PC_VIZIER}", version=description, rows=len(table))
    return table


# --------------------------------------------------------------------------
# 5. Mamajek's dwarf sequence
# --------------------------------------------------------------------------

def fetch_dwarf_sequence(*, refresh: bool = False) -> str:
    """Download "A Modern Mean Dwarf Stellar Color and Effective Temperature Sequence"."""
    if C.DWARF_SEQUENCE_PATH.exists() and not refresh:
        return C.DWARF_SEQUENCE_PATH.read_text()
    print("  Mamajek dwarf colour-temperature sequence")
    with urllib.request.urlopen(DWARF_SEQUENCE_URL, timeout=120) as r:
        text = r.read().decode("utf-8", errors="replace")
    C.DWARF_SEQUENCE_PATH.write_text(text)
    version = re.search(r"# Version (\S+)", text)
    _record_provenance("dwarf_sequence", source="Mamajek, 'A Modern Mean Dwarf Stellar Color "
                       "and Effective Temperature Sequence' (see Pecaut & Mamajek 2013, ApJS 208, 9)",
                       url=DWARF_SEQUENCE_URL, version=version.group(1) if version else None)
    return text


# --------------------------------------------------------------------------
# 6. SIMBAD: names and famous stars
# --------------------------------------------------------------------------

def _simbad():
    from astroquery.simbad import Simbad

    s = Simbad()
    s.add_votable_fields("ids", "plx_value", "plx_bibcode", "V", "sp_type", "otype")
    return s


def _parse_ids(ids: str) -> dict:
    """Pull the Gaia DR3 id, HIP number and common names out of SIMBAD's id list."""
    out = {"gaia_dr3": None, "hip": None, "names": []}
    for ident in str(ids).split("|"):
        ident = ident.strip()
        if ident.startswith("Gaia DR3 "):
            out["gaia_dr3"] = int(ident.split()[-1])
        elif ident.startswith("HIP ") and ident[4:].strip().isdigit():
            out["hip"] = int(ident[4:])
        elif ident.startswith("NAME "):
            out["names"].append(ident[5:])
    return out


def _float_or_none(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if not np.isfinite(f) else f


def fetch_famous(*, refresh: bool = False) -> dict:
    """Resolve every star in config.FAMOUS_STARS / BEYOND_STARS, plus the clusters."""
    if C.FAMOUS_PATH.exists() and not refresh:
        return json.loads(C.FAMOUS_PATH.read_text())
    print("  SIMBAD: resolving famous stars and reference clusters")
    stars = C.FAMOUS_STARS + C.BEYOND_STARS
    names = [s.simbad for s in stars] + list(C.REFERENCE_CLUSTERS.values())
    result = _simbad().query_objects(names)
    out = {}
    for row in result:
        # object_number_id is the 1-based position in the request. (The echoed
        # user_specified_id is space-padded, so it is not a safe key.)
        ident = names[int(row["object_number_id"]) - 1]
        parsed = _parse_ids(row["ids"])
        out[ident] = {
            "main_id": str(row["main_id"]),
            "ra": _float_or_none(row["ra"]), "dec": _float_or_none(row["dec"]),
            "plx_mas": _float_or_none(row["plx_value"]),
            "plx_bibcode": str(row["plx_bibcode"]) or None,
            "vmag": _float_or_none(row["V"]),
            "sp_type": str(row["sp_type"]) or None,
            "otype": str(row["otype"]),
            **parsed,
        }
    missing = [n for n in names if n not in out]
    if missing:
        raise RuntimeError(f"SIMBAD could not resolve: {missing}")
    C.FAMOUS_PATH.write_text(json.dumps(out, indent=2))
    _record_provenance("simbad_famous", source="SIMBAD (CDS, Strasbourg), query_objects",
                       objects=names)
    return out


def fetch_names(gaia: Table, hip: Table, *, refresh: bool = False) -> Table:
    """SIMBAD main identifiers and common names for every star within 20 pc.

    Uploads the Gaia DR3 and HIP identifiers of all candidates within the named
    radius and joins them against SIMBAD's identifier table in one query each.
    SIMBAD spectral types come along too, as an independent check on our own
    colour-based classification.
    """
    if C.NAMES_PATH.exists() and not refresh:
        return Table.read(C.NAMES_PATH)
    print(f"  SIMBAD: names for every star within {C.NAMED_RADIUS_PC:g} pc")
    s = _simbad()
    min_plx = 1000.0 / C.NAMED_RADIUS_PC
    ids = [f"Gaia DR3 {sid}" for sid in gaia["source_id"][gaia["parallax"] > min_plx]]
    ids += [f"HIP {h}" for h in hip["hip"][hip["plx"] > min_plx]]
    upload = Table({"ident": ids})

    basic = s.query_tap(
        "SELECT u.ident, b.main_id, b.otype, b.sp_type, b.plx_value "
        "FROM TAP_UPLOAD.up AS u JOIN ident AS i ON i.id = u.ident "
        "JOIN basic AS b ON b.oid = i.oidref", up=upload, maxrec=100000)
    common = s.query_tap(
        "SELECT u.ident, n.id AS name FROM TAP_UPLOAD.up AS u "
        "JOIN ident AS i ON i.id = u.ident JOIN ident AS n ON n.oidref = i.oidref "
        "WHERE n.id LIKE 'NAME %'", up=upload, maxrec=100000)

    names: dict[str, list[str]] = {}
    for row in common:
        names.setdefault(str(row["ident"]), []).append(str(row["name"])[5:])
    table = Table({
        "ident": [str(x) for x in basic["ident"]],
        "main_id": [str(x) for x in basic["main_id"]],
        "otype": [str(x) for x in basic["otype"]],
        "sp_type": [str(x) if x is not np.ma.masked else "" for x in basic["sp_type"]],
        "common_names": ["|".join(names.get(str(x), [])) for x in basic["ident"]],
    })
    table.write(C.NAMES_PATH, overwrite=True)
    print(f"  {len(table):,} of {len(ids):,} identifiers matched in SIMBAD")
    _record_provenance("simbad_names", source="SIMBAD (CDS) TAP, identifier upload join",
                       uploaded=len(ids), matched=len(table))
    return table


# --------------------------------------------------------------------------
# 7. Known planet hosts
# --------------------------------------------------------------------------

def fetch_exoplanet_hosts(*, refresh: bool = False) -> Table:
    """Every confirmed planet whose host lies within the map, one row per planet."""
    if C.EXOPLANET_HOSTS_PATH.exists() and not refresh:
        return Table.read(C.EXOPLANET_HOSTS_PATH)
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

    print(f"  NASA Exoplanet Archive: planets with hosts within {C.MAX_DISTANCE_PC:g} pc")
    table = NasaExoplanetArchive.query_criteria(
        table="pscomppars", select=",".join(EXOPLANET_COLUMNS),
        where=f"sy_dist < {C.MAX_DISTANCE_PC:g}")
    plain = Table()
    for name in EXOPLANET_COLUMNS:  # astroquery also attaches a SkyCoord; skip it
        col = table[name]
        values = np.ma.asarray(getattr(col, "value", col))  # drop units FITS cannot store
        if values.dtype.kind in "fiu":
            plain[name] = np.ma.filled(values.astype(float), np.nan)
        else:
            plain[name] = [str(v) if v is not np.ma.masked else "" for v in values]
    plain = _fix_strings(plain)
    plain.write(C.EXOPLANET_HOSTS_PATH, overwrite=True)
    print(f"  {len(plain):,} planets around {len(set(plain['hostname'])):,} host stars")
    _record_provenance("exoplanet_hosts", source="NASA Exoplanet Archive, pscomppars",
                       where=f"sy_dist < {C.MAX_DISTANCE_PC:g}", rows=len(plain))
    return plain


# --------------------------------------------------------------------------

def fetch_all(*, refresh: bool = False) -> None:
    gaia = fetch_gaia(refresh=refresh)
    fetch_gcns(refresh=refresh)
    hip = fetch_hipparcos(refresh=refresh)
    fetch_naked_eye(refresh=refresh)
    fetch_census_10pc(refresh=refresh)
    fetch_dwarf_sequence(refresh=refresh)
    fetch_famous(refresh=refresh)
    fetch_names(gaia, hip, refresh=refresh)
    fetch_exoplanet_hosts(refresh=refresh)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download everything, ignoring the cache")
    args = parser.parse_args(argv)
    print("Fetching data for the solar-neighbourhood map")
    fetch_all(refresh=args.refresh)
    print("\nDone. Cache is warm; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
