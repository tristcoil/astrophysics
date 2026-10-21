"""Download and cache everything the analysis needs.

Four things are fetched, all free and public:

1. **Every Gaia DR3 source with parallax > 10 mas** from the ESA Gaia archive --
   the candidates for "stars within 100 parsecs". About 575,000 rows. No
   quality cuts are applied on the server: we want to *see* what each cut
   removes, so the cuts happen locally in ``analysis.py``.

2. **Star counts per sky pixel** for the whole Gaia catalogue, counted by the
   archive from a random 2% subsample (see ``config.SKY_QUERY`` for why).

3. **Eric Mamajek's mean dwarf sequence**, a published table linking spectral
   type, surface temperature and Gaia colour for main-sequence stars. It lets
   the diagram's colour axis be read as a temperature.

4. **Gaia identifiers for a handful of famous stars**, resolved through the
   SIMBAD database at CDS, so the diagram can say "this dot is Proxima
   Centauri".

Everything is cached under ``data/`` with its retrieval date and the exact query
that produced it. Nothing is re-downloaded once present, so ``run.py`` works
offline after one run of this script::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.request
import warnings
from datetime import datetime, timezone

import numpy as np
from astropy.table import Table

from config import (
    BRIGHT_STARS_ABSENT,
    CANDIDATE_QUERY,
    CANDIDATES_PATH,
    DWARF_SEQUENCE_PATH,
    NAMED_STARS,
    NAMED_STARS_PATH,
    SKY_COUNTS_PATH,
    SKY_HEALPIX_LEVEL,
    SKY_QUERY,
    SKY_SAMPLE_SIZE,
)

warnings.filterwarnings("ignore", category=UserWarning)

DWARF_SEQUENCE_URL = (
    "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt"
)

#: TAP services holding the full Gaia DR3 gaia_source table. ESA's archive is
#: the canonical source; Gaia@AIP (Leibniz-Institut fuer Astrophysik Potsdam) is
#: an official mirror of the same release, and can be far faster when ESA's
#: archive is busy. Both return 574,531 rows for the candidate query.
TAP_SERVICES = {
    "esa": ("ESA Gaia archive", "https://gea.esac.esa.int/tap-server/tap"),
    "aip": ("Gaia@AIP mirror", "https://gaia.aip.de/tap"),
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _run_async(query: str, *, label: str, service: str = "esa"):
    """Run an ADQL query as an asynchronous TAP job and return the table.

    Asynchronous jobs matter: the Gaia archive silently truncates synchronous
    queries at 2,000 rows, which would be an easy way to lose most of a sample
    without any error.
    """
    import pyvo

    name, url = TAP_SERVICES[service]
    print(f"  submitting {label} to the {name} (asynchronous job)...")
    t0 = time.time()
    job = pyvo.dal.TAPService(url).submit_job(query, maxrec=5_000_000)
    job.run()
    # Logged before waiting: if this process dies after the service finishes,
    # the result can still be collected by the job's URL.
    print(f"  job {job.url}; waiting for it to run...", flush=True)
    job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=4 * 3600)
    job.raise_if_error()
    table = job.fetch_result().to_table()
    print(f"  {len(table):,} rows in {time.time() - t0:.0f} s")
    return table


def _count(query_where: str, *, service: str) -> int:
    """Row count the server expects, to prove a download is complete."""
    import pyvo

    q = f"SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE {query_where}"
    return int(pyvo.dal.TAPService(TAP_SERVICES[service][1]).run_sync(q).to_table()["n"][0])


def _unmask(table: Table) -> Table:
    """A plain table with NaN wherever the archive returned a null.

    The archive returns *every* column as a masked column, including integer
    ones such as source_id, which cannot hold NaN -- so a blanket
    ``table.filled(np.nan)`` fails. Nulls are expected only in float columns (a
    source with no BP/RP photometry, no radial velocity); a null in an integer
    column would be a surprise worth stopping for.
    """
    from astropy.table import Column

    columns = []
    for name in table.colnames:
        col = table[name]
        mask = np.asarray(getattr(col, "mask", False), dtype=bool)
        data = np.array(col, copy=True) if not hasattr(col, "mask") else np.array(col.data)
        if col.dtype.kind == "f":
            data = np.where(mask, np.nan, data)
        elif mask.any():
            raise RuntimeError(f"unexpected nulls in integer column {name!r}")
        columns.append(Column(data, name=name, unit=getattr(col, "unit", None),
                              description=getattr(col, "description", None)))
    return Table(columns)


def _write_with_provenance(table: Table, path, *, query: str, source: str,
                           retrieved: str | None = None, **extra) -> None:
    """Write a FITS table whose header records where the data came from."""
    table.meta.clear()
    table.meta["SOURCE"] = source
    table.meta["RETRIEVD"] = retrieved or _now()
    # FITS header values must be printable ASCII on one line.
    table.meta["QUERY"] = " ".join(query.split())
    for key, value in extra.items():
        table.meta[key.upper()[:8]] = value
    table.write(path, format="fits", overwrite=True)
    size_mb = path.stat().st_size / 1e6
    print(f"  wrote {path.relative_to(path.parent.parent)} ({size_mb:.1f} MB)")


# --------------------------------------------------------------------------
# 1. The 100 pc candidates
# --------------------------------------------------------------------------

def fetch_candidates(*, refresh: bool = False, service: str = "aip") -> None:
    if CANDIDATES_PATH.exists() and not refresh:
        print(f"  cached: {CANDIDATES_PATH.name}")
        return
    expected = _count(CANDIDATE_QUERY.split("WHERE", 1)[1].strip(), service=service)
    table = _unmask(_run_async(CANDIDATE_QUERY, label="the parallax > 10 mas query",
                               service=service))
    # A TAP service may cap a result at its row limit and say so only in a
    # warning. Comparing with the server's own count makes truncation fatal.
    if len(table) != expected:
        raise RuntimeError(f"downloaded {len(table):,} rows but the server counts "
                           f"{expected:,}: the result was truncated")
    name, url = TAP_SERVICES[service]
    _write_with_provenance(
        table, CANDIDATES_PATH, query=CANDIDATE_QUERY,
        source=f"Gaia DR3 gaiadr3.gaia_source via the {name} ({url})",
    )


def load_candidates() -> Table:
    if not CANDIDATES_PATH.exists():
        raise FileNotFoundError(
            f"{CANDIDATES_PATH} not found -- run `python fetch_data.py` first"
        )
    # mask_invalid=False: keep NaN as NaN rather than hiding it behind a mask.
    return Table.read(CANDIDATES_PATH, format="fits", mask_invalid=False)


# --------------------------------------------------------------------------
# 2. The all-sky star counts
# --------------------------------------------------------------------------

def _catalogue_size() -> int:
    """Number of sources in Gaia DR3, from the top of its ``random_index``.

    ``random_index`` is a permutation of 0..N-1, so its maximum plus one is the
    catalogue size. Asking for the maximum uses the column's index and returns
    in seconds, whereas COUNT(*) would scan all 1.8 billion rows.
    """
    from astroquery.gaia import Gaia

    job = Gaia.launch_job(
        "SELECT TOP 1 random_index FROM gaiadr3.gaia_source_lite "
        "ORDER BY random_index DESC"
    )
    return int(job.get_results()["random_index"][0]) + 1


def fetch_sky_counts(*, refresh: bool = False) -> None:
    if SKY_COUNTS_PATH.exists() and not refresh:
        print(f"  cached: {SKY_COUNTS_PATH.name}")
        return
    n_total = _catalogue_size()
    print(f"  Gaia DR3 holds {n_total:,} sources; counting a random "
          f"{SKY_SAMPLE_SIZE:,} ({SKY_SAMPLE_SIZE / n_total:.2%}) by sky pixel. "
          f"This takes several minutes.")
    table = _run_async(SKY_QUERY, label="the sky-count query")
    write_sky_counts(table, n_total=n_total)


def write_sky_counts(table: Table, *, n_total: int, retrieved: str | None = None):
    """Store the per-pixel counts with everything needed to scale them up."""
    table = Table({"hpx": np.asarray(table[table.colnames[0]], dtype=np.int64),
                   "n": np.asarray(table["n"], dtype=np.int64)})
    counted = int(table["n"].sum())
    if counted != SKY_SAMPLE_SIZE:
        raise RuntimeError(
            f"counts sum to {counted:,}, expected {SKY_SAMPLE_SIZE:,}: the "
            f"subsample is incomplete"
        )
    _write_with_provenance(
        table, SKY_COUNTS_PATH, query=SKY_QUERY, retrieved=retrieved,
        source="ESA Gaia archive, gaiadr3.gaia_source_lite",
        level=SKY_HEALPIX_LEVEL, nsample=SKY_SAMPLE_SIZE, ntotal=n_total,
    )


def load_sky_counts() -> Table:
    if not SKY_COUNTS_PATH.exists():
        raise FileNotFoundError(
            f"{SKY_COUNTS_PATH} not found -- run `python fetch_data.py` first"
        )
    return Table.read(SKY_COUNTS_PATH, format="fits")


# --------------------------------------------------------------------------
# 3. Mamajek's dwarf sequence
# --------------------------------------------------------------------------

def fetch_dwarf_sequence(*, refresh: bool = False) -> None:
    if DWARF_SEQUENCE_PATH.exists() and not refresh:
        print(f"  cached: {DWARF_SEQUENCE_PATH.name}")
        return
    print(f"  downloading {DWARF_SEQUENCE_URL}")
    request = urllib.request.Request(DWARF_SEQUENCE_URL,
                                     headers={"User-Agent": "astrophysics-atlas"})
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8")
    header = f"# Retrieved {_now()} from {DWARF_SEQUENCE_URL}\n"
    DWARF_SEQUENCE_PATH.write_text(header + text)
    print(f"  wrote {DWARF_SEQUENCE_PATH.name}")


def load_dwarf_sequence() -> dict[str, np.ndarray]:
    """Parse the main table of Mamajek's file into arrays keyed by column name.

    The file is whitespace-aligned text with ``...`` for "no value". Its first
    block, introduced by a ``#SpT`` header line, is the dwarf sequence proper;
    the block ends at the next comment line.
    """
    if not DWARF_SEQUENCE_PATH.exists():
        raise FileNotFoundError(
            f"{DWARF_SEQUENCE_PATH} not found -- run `python fetch_data.py` first"
        )
    lines = DWARF_SEQUENCE_PATH.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("#SpT"))
    names = lines[start].lstrip("#").split()
    rows = []
    for line in lines[start + 1:]:
        if line.startswith("#") or not line.strip():
            break
        rows.append(line.split())

    def number(v: str) -> float:
        # A trailing colon marks a value the table flags as uncertain.
        try:
            return float(v.rstrip(":"))
        except ValueError:
            return np.nan

    columns: dict[str, np.ndarray] = {}
    for j, name in enumerate(names):
        if name in columns:  # the file repeats SpT as its last column
            continue
        values = [r[j] if j < len(r) else "..." for r in rows]
        if name == "SpT":
            columns[name] = np.array(values)
        else:
            columns[name] = np.array([number(v) for v in values])
    return columns


# --------------------------------------------------------------------------
# 4. Famous stars
# --------------------------------------------------------------------------

def fetch_named_stars(*, refresh: bool = False) -> None:
    """Resolve each named star's Gaia DR3 source_id through SIMBAD.

    Also checks, with a one-arcminute cone search, that each "absent" bright
    star really has no Gaia DR3 counterpart, recording the brightest source that
    *is* there so the claim can be audited.
    """
    if NAMED_STARS_PATH.exists() and not refresh:
        print(f"  cached: {NAMED_STARS_PATH.name}")
        return
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia
    from astroquery.simbad import Simbad

    out = {"_retrieved_utc": _now(), "_source": "SIMBAD (CDS) identifiers",
           "stars": {}, "absent_bright_stars": {}}

    for display, simbad_id in NAMED_STARS.items():
        ids = [str(i) for i in Simbad.query_objectids(simbad_id)["id"]]
        gaia = [i for i in ids if i.startswith("Gaia DR3 ")]
        source_id = int(gaia[0].split()[-1]) if gaia else None
        out["stars"][display] = {"simbad": simbad_id, "source_id": source_id}
        print(f"  {display:<20} {source_id}")

    for name in BRIGHT_STARS_ABSENT:
        c = SkyCoord.from_name(name)
        query = (
            "SELECT TOP 1 source_id, phot_g_mean_mag, parallax, "
            f"DISTANCE(POINT({c.ra.deg}, {c.dec.deg}), POINT(ra, dec)) * 3600 "
            "AS sep_arcsec FROM gaiadr3.gaia_source "
            f"WHERE 1 = CONTAINS(POINT(ra, dec), CIRCLE({c.ra.deg}, {c.dec.deg}, "
            f"{1 / 60})) ORDER BY phot_g_mean_mag ASC"
        )
        r = Gaia.launch_job(query).get_results()
        brightest = None
        if len(r):
            brightest = {"source_id": int(r["source_id"][0]),
                         "g_mag": float(r["phot_g_mean_mag"][0]),
                         "sep_arcsec": float(r["sep_arcsec"][0])}
        out["absent_bright_stars"][name] = {
            "ra": c.ra.deg, "dec": c.dec.deg,
            "brightest_dr3_source_within_1arcmin": brightest,
        }
        g = f"G = {brightest['g_mag']:.2f} at {brightest['sep_arcsec']:.0f}\"" \
            if brightest else "nothing"
        print(f"  {name:<20} absent from DR3; brightest source within 1': {g}")

    NAMED_STARS_PATH.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"  wrote {NAMED_STARS_PATH.name}")


def load_named_stars() -> dict:
    if not NAMED_STARS_PATH.exists():
        raise FileNotFoundError(
            f"{NAMED_STARS_PATH} not found -- run `python fetch_data.py` first"
        )
    return json.loads(NAMED_STARS_PATH.read_text())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download everything, ignoring the cache")
    parser.add_argument("--skip-sky", action="store_true",
                        help="skip the slow all-sky star-count query")
    parser.add_argument("--mirror", choices=sorted(TAP_SERVICES), default="aip",
                        help="where to fetch the 100 pc sample from (default: aip, "
                             "which answered in 72 s when ESA's archive took over an hour)")
    args = parser.parse_args(argv)

    print("1. Gaia DR3 sources with parallax > 10 mas")
    fetch_candidates(refresh=args.refresh, service=args.mirror)
    if not args.skip_sky:
        print("\n2. All-sky star counts")
        fetch_sky_counts(refresh=args.refresh)
    print("\n3. Mamajek's mean dwarf sequence")
    fetch_dwarf_sequence(refresh=args.refresh)
    print("\n4. Famous stars")
    fetch_named_stars(refresh=args.refresh)

    print("\nDone. Cache is warm; run.py can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
