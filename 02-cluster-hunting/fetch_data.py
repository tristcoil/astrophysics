"""Download and cache everything the analysis needs.

Three kinds of data are fetched, all free and public:

1. **Gaia DR3 astrometry and photometry** for every star brighter than G = 18
   in a cone around each cluster: positions, parallaxes, proper motions, their
   uncertainties and correlations, G/BP/RP photometry, and radial velocities
   where Gaia measured them. This is the only data the analysis uses to find
   and characterise the clusters.

   The table is read from the CDS VizieR mirror of the Gaia DR3 catalogue
   (I/355/gaiadr3) rather than from ESA's own Gaia archive. The values are
   identical -- VizieR serves a copy of the official release -- but at the
   time of writing ESA's archive was timing out while it is rebuilt for DR4,
   whereas VizieR returned the Pleiades field in under 20 seconds. Columns are
   renamed to the official Gaia archive names on the way in, so the rest of the
   code reads exactly like any other Gaia analysis. ``--service esa`` switches
   to ESA's archive when it is healthy.

2. **Published cluster catalogues**, used only at the end to grade the result:
   Hunt & Reffert (2023), the most complete Gaia DR3 census of star clusters,
   including their member lists keyed by Gaia source ID; and the Gaia
   Collaboration's (2018) DR2 cluster paper, which lists each nearby cluster's
   3D space velocity.

3. **Names of famous stars** (the Seven Sisters themselves) resolved through
   SIMBAD to Gaia DR3 source IDs, purely so the figures can label them.

Everything is cached under ``data/`` with its retrieval date and the exact query
used. Once the cache is warm the analysis runs with no network at all.

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from config import CLUSTERS, GAIA_DIR, LITERATURE_DIR, Cluster, cluster_by_key

warnings.filterwarnings("ignore")

VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
VIZIER_GAIA_TABLE = '"I/355/gaiadr3"'

#: VizieR column name -> official Gaia DR3 archive column name.
GAIA_COLUMNS = {
    "Source": "source_id",
    "RA_ICRS": "ra",
    "DE_ICRS": "dec",
    "GLON": "l",
    "GLAT": "b",
    "ELAT": "ecl_lat",
    "Plx": "parallax",
    "e_Plx": "parallax_error",
    "pmRA": "pmra",
    "e_pmRA": "pmra_error",
    "pmDE": "pmdec",
    "e_pmDE": "pmdec_error",
    "PlxpmRAcor": "parallax_pmra_corr",
    "PlxpmDEcor": "parallax_pmdec_corr",
    "pmRApmDEcor": "pmra_pmdec_corr",
    "RUWE": "ruwe",
    "Solved": "astrometric_params_solved",
    "nueff": "nu_eff_used_in_astrometry",
    "pscol": "pseudocolour",
    "Gmag": "phot_g_mean_mag",
    "BPmag": "phot_bp_mean_mag",
    "RPmag": "phot_rp_mean_mag",
    "BP-RP": "bp_rp",
    "E(BP/RP)": "phot_bp_rp_excess_factor",
    "RV": "radial_velocity",
    "e_RV": "radial_velocity_error",
    "GRVSmag": "grvs_mag",
    "RVtempTeff": "rv_template_teff",
}

HR23_CLUSTER_COLUMNS = [
    "Name", "Type", "N", "Nt", "RA_ICRS", "DE_ICRS", "r50", "rc", "rt", "rtot",
    "r50pc", "rtpc", "pmRA", "s_pmRA", "e_pmRA", "pmDE", "s_pmDE", "e_pmDE",
    "Plx", "s_Plx", "e_Plx", "dist16", "dist50", "dist84", "RV", "s_RV", "e_RV",
    "n_RV", "logAge16", "logAge50", "logAge84", "AV16", "AV50", "AV84", "MOD50",
]
HR23_MEMBER_COLUMNS = ["GaiaDR3", "Prob", "inrt", "Gmag", "RA_ICRS", "DE_ICRS"]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _tap_sync(query: str, *, attempts: int = 3, maxrec: int = 3_000_000):
    """Run a synchronous VizieR TAP query, retrying transient failures."""
    import pyvo

    service = pyvo.dal.TAPService(VIZIER_TAP)
    for attempt in range(1, attempts + 1):
        try:
            return service.run_sync(query, maxrec=maxrec, timeout=900).to_table()
        except Exception as exc:  # network hiccups, 5xx from the service
            if attempt == attempts:
                raise
            wait = 10 * attempt
            print(f"    query failed ({type(exc).__name__}); retrying in {wait}s")
            time.sleep(wait)


def _plain(table) -> pd.DataFrame:
    """Astropy table -> DataFrame with NaN for every missing float value.

    Gaia leaves many fields legitimately empty -- no radial velocity for faint
    stars, no pseudocolour for 5-parameter solutions -- and VizieR returns
    those as masked entries. Converting once, here, means no downstream code
    has to think about masks.
    """
    out = {}
    for name in table.colnames:
        col = table[name]
        data = np.ma.getdata(col)
        mask = np.ma.getmaskarray(col)
        if data.dtype.kind == "f":
            values = np.asarray(data, dtype=float).copy()
            values[mask] = np.nan
        elif data.dtype.kind in "iu":
            # FITS stores integers big-endian; pandas indexing needs native order.
            values = np.asarray(data).astype(np.int64)
            if mask.any():
                values = np.where(mask, -1, values)
        else:
            values = np.asarray(data).astype(str)
        out[name] = values
    return pd.DataFrame(out)


# --------------------------------------------------------------------------
# Gaia DR3
# --------------------------------------------------------------------------

def gaia_paths(cluster: Cluster):
    stem = GAIA_DIR / f"{cluster.key}_gaia_dr3"
    return stem.with_suffix(".fits"), stem.with_suffix(".json")


def gaia_query(cluster: Cluster, service: str = "vizier") -> str:
    """The ADQL for one cluster's cone.

    Only two cuts are made at the archive: the magnitude limit, and requiring a
    parallax, which keeps the 5- and 6-parameter astrometric solutions and
    drops the 2-parameter ones that have a position but no measured motion.
    Every other quality decision is made locally, where it can be inspected and
    changed without re-downloading.
    """
    if service == "vizier":
        cols = ", ".join(f'"{c}"' for c in GAIA_COLUMNS)
        return (
            f"SELECT {cols} FROM {VIZIER_GAIA_TABLE}\n"
            f"WHERE 1 = CONTAINS(POINT('ICRS', \"RA_ICRS\", \"DE_ICRS\"),\n"
            f"                   CIRCLE('ICRS', {cluster.ra_deg}, {cluster.dec_deg}, "
            f"{cluster.radius_deg}))\n"
            f"  AND \"Gmag\" < {cluster.g_max} AND \"Plx\" IS NOT NULL"
        )
    if service == "esa":
        cols = ", ".join(GAIA_COLUMNS.values())
        return (
            f"SELECT {cols} FROM gaiadr3.gaia_source\n"
            f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),\n"
            f"                   CIRCLE('ICRS', {cluster.ra_deg}, {cluster.dec_deg}, "
            f"{cluster.radius_deg}))\n"
            f"  AND phot_g_mean_mag < {cluster.g_max} AND parallax IS NOT NULL"
        )
    raise ValueError(f"unknown service {service!r}")


def fetch_gaia(cluster: Cluster, *, service: str = "vizier", refresh: bool = False,
               quiet: bool = False) -> pd.DataFrame:
    """Download (or read from cache) the Gaia DR3 cone for one cluster."""
    fits_path, meta_path = gaia_paths(cluster)
    if fits_path.exists() and not refresh:
        return load_gaia(cluster)

    query = gaia_query(cluster, service)
    if not quiet:
        print(f"  querying Gaia DR3 via {service}: {cluster.radius_deg} deg around "
              f"({cluster.ra_deg}, {cluster.dec_deg}), G < {cluster.g_max}")
    t0 = time.time()
    if service == "vizier":
        table = _tap_sync(query)
        for old, new in GAIA_COLUMNS.items():
            table.rename_column(old, new)
    else:
        from astroquery.gaia import Gaia

        table = Gaia.launch_job_async(query).get_results()
        for name in list(table.colnames):
            table.rename_column(name, name.lower())

    table.write(fits_path, overwrite=True)
    meta_path.write_text(json.dumps({
        "_retrieved_utc": _now(),
        "_service": VIZIER_TAP if service == "vizier" else "ESA Gaia archive",
        "_catalogue": "Gaia DR3 (Gaia Collaboration, Vallenari et al. 2023)",
        "query": query,
        "rows": len(table),
    }, indent=2))
    if not quiet:
        size_mb = fits_path.stat().st_size / 1e6
        print(f"  {len(table):,} stars in {time.time() - t0:.0f}s -> "
              f"{fits_path.name} ({size_mb:.1f} MB)")
    return load_gaia(cluster)


def load_gaia(cluster: Cluster) -> pd.DataFrame:
    """Read the cached cone, erroring clearly if it has not been fetched."""
    from astropy.table import Table

    fits_path, _ = gaia_paths(cluster)
    if not fits_path.exists():
        raise FileNotFoundError(f"{fits_path} not found -- run `python fetch_data.py` first")
    df = _plain(Table.read(fits_path))
    df["source_id"] = df["source_id"].astype(np.int64)
    return df


def gaia_provenance(cluster: Cluster) -> dict:
    _, meta_path = gaia_paths(cluster)
    return json.loads(meta_path.read_text()) if meta_path.exists() else {}


# --------------------------------------------------------------------------
# Published catalogues
# --------------------------------------------------------------------------

LITERATURE_PATH = LITERATURE_DIR / "clusters.json"


def _member_path(cluster: Cluster):
    return LITERATURE_DIR / f"hr23_members_{cluster.key}.csv"


def _row_dict(table, row) -> dict:
    out = {}
    for name in table.colnames:
        value = row[name]
        if np.ma.is_masked(value):
            out[name] = None
            continue
        if hasattr(value, "item"):
            value = value.item()
        if isinstance(value, bytes):
            value = value.decode()
        if isinstance(value, float) and value != value:
            value = None
        out[name] = value
    return out


def fetch_literature(clusters: list[Cluster], *, quiet: bool = False) -> dict:
    """Published parameters and member lists, cached with the retrieval date."""
    out = {
        "_retrieved_utc": _now(),
        "_sources": {
            "hr23": "Hunt & Reffert (2023), A&A 673, A114 -- VizieR J/A+A/673/A114",
            "gc2018": (
                "Gaia Collaboration, Babusiaux et al. (2018), A&A 616, A10 -- "
                "VizieR J/A+A/616/A10 (Gaia DR2; U,V,W are ICRS Cartesian components)"
            ),
        },
        "clusters": {},
    }

    names = ", ".join(f"'{c.catalogue_name}'" for c in clusters)
    cols = ", ".join(f'"{c}"' for c in HR23_CLUSTER_COLUMNS)
    hr23 = _tap_sync(f'SELECT {cols} FROM "J/A+A/673/A114/clusters" WHERE "Name" IN ({names})')
    hr23_rows = {str(r["Name"]): _row_dict(hr23, r) for r in hr23}

    gc_near = _tap_sync('SELECT * FROM "J/A+A/616/A10/tablea3"')
    gc_far = _tap_sync('SELECT * FROM "J/A+A/616/A10/tablea4"')
    gc_rows = {str(r["Cluster"]).strip(): _row_dict(gc_near, r) for r in gc_near}
    gc_rows.update({str(r["Cluster"]).strip(): _row_dict(gc_far, r) for r in gc_far})

    for cluster in clusters:
        if not quiet:
            print(f"  {cluster.name}: Hunt & Reffert (2023) + Gaia Collaboration (2018)")
        entry = {
            "hr23": hr23_rows.get(cluster.catalogue_name),
            "gc2018": gc_rows.get(cluster.gc2018_name),
        }
        if entry["hr23"] is None:
            print(f"    warning: {cluster.catalogue_name} not found in HR23")
        if entry["gc2018"] is None:
            print(f"    warning: {cluster.gc2018_name} not found in GC2018")
        out["clusters"][cluster.key] = entry

        mcols = ", ".join(f'"{c}"' for c in HR23_MEMBER_COLUMNS)
        members = _tap_sync(
            f'SELECT {mcols} FROM "J/A+A/673/A114/members" '
            f"WHERE \"Name\" = '{cluster.catalogue_name}'"
        )
        df = _plain(members).rename(columns={"GaiaDR3": "source_id"})
        df["source_id"] = df["source_id"].astype(np.int64)
        df.to_csv(_member_path(cluster), index=False)
        if not quiet:
            print(f"    {len(df):,} published members cached")

        famous = resolve_famous_stars(cluster)
        out["clusters"][cluster.key]["famous_stars"] = famous
        _fetch_named_rows(cluster, famous)

    LITERATURE_PATH.write_text(json.dumps(out, indent=2))
    return out


def _named_path(cluster: Cluster):
    return LITERATURE_DIR / f"named_stars_{cluster.key}.csv"


def _fetch_named_rows(cluster: Cluster, famous: dict) -> None:
    """Gaia rows for the named stars, *including* those with no parallax.

    The main cone requires a parallax, which silently drops stars Gaia could
    only position -- among them the brightest of the Pleiades, Alcyone. Fetching
    them by ID lets the figures show them and say why they are not selected.
    """
    ids = [str(i) for i in famous.values() if i]
    if not ids:
        return
    cols = ", ".join(f'"{c}"' for c in GAIA_COLUMNS)
    table = _tap_sync(f'SELECT {cols} FROM {VIZIER_GAIA_TABLE} WHERE "Source" IN ({", ".join(ids)})')
    for old, new in GAIA_COLUMNS.items():
        table.rename_column(old, new)
    df = _plain(table)
    names = {v: k for k, v in famous.items() if v}
    df["name"] = [names.get(int(s), "") for s in df["source_id"]]
    df.to_csv(_named_path(cluster), index=False)


def load_named_stars(cluster: Cluster) -> pd.DataFrame | None:
    path = _named_path(cluster)
    return pd.read_csv(path) if path.exists() else None


def resolve_famous_stars(cluster: Cluster) -> dict:
    """Named stars -> Gaia DR3 source IDs, via SIMBAD. Labels only."""
    if not cluster.famous_stars:
        return {}
    from astroquery.simbad import Simbad

    found = {}
    for name in cluster.famous_stars:
        try:
            ids = Simbad.query_objectids(name)
            col = [c for c in ids.colnames if c.lower() in ("id", "main_id")][0]
            gaia = [str(x) for x in ids[col] if str(x).startswith("Gaia DR3 ")]
            found[name] = int(gaia[0].split()[-1]) if gaia else None
        except Exception as exc:
            print(f"    could not resolve {name}: {exc}")
            found[name] = None
    return found


def load_literature() -> dict:
    if not LITERATURE_PATH.exists():
        raise FileNotFoundError(f"{LITERATURE_PATH} not found -- run `python fetch_data.py` first")
    return json.loads(LITERATURE_PATH.read_text())


def load_published_members(cluster: Cluster) -> pd.DataFrame | None:
    path = _member_path(cluster)
    return pd.read_csv(path) if path.exists() else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--cluster", action="append",
                        help="fetch only this cluster (repeatable); default is all")
    parser.add_argument("--service", choices=("vizier", "esa"), default="vizier",
                        help="where to read Gaia DR3 from (default: CDS VizieR)")
    parser.add_argument("--refresh", action="store_true",
                        help="re-download even if cached")
    parser.add_argument("--skip-literature", action="store_true",
                        help="skip the published catalogues")
    args = parser.parse_args(argv)

    clusters = [cluster_by_key(k) for k in args.cluster] if args.cluster else CLUSTERS

    print("Fetching Gaia DR3 cones")
    for cluster in clusters:
        print(f"\n{cluster.name} -- {cluster.difficulty} case")
        df = fetch_gaia(cluster, service=args.service, refresh=args.refresh)
        print(f"  {len(df):,} stars cached")

    if not args.skip_literature:
        print("\nFetching published catalogues (used only for grading)")
        fetch_literature(CLUSTERS)

    print("\nDone. Cache is warm; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
