"""Download and cache everything the analysis needs.

Four kinds of data, all public:

1. **THINGS** (Walter et al. 2008): VLA 21-cm data cubes of four galaxies --
   NGC 3198, NGC 2403, NGC 2903 and DDO 154 -- plus the survey's own moment
   maps. A cube is the raw material of this project: an image of the sky at
   each of 57-87 radio frequencies, so every pixel holds a spectrum of the
   hydrogen line whose Doppler shift says how fast the gas there moves toward
   or away from us. About 2 GB, mostly NGC 2403's cube.

2. **Published results**, used only to grade ours: the LaTeX and figure
   sources of de Blok et al. (2008), Trachternach et al. (2008) and Oh et al.
   (2011) from arXiv (their rotation curves exist only as figures, which are
   vector PostScript, so the plotted points can be read back exactly), the
   THINGS galaxy table and global HI spectra from VizieR, and SPARC's rotation
   curves of the same galaxies from older telescopes (project D4).

3. **MaNGA DR17 catalogues**: the DAP summary table (one row per galaxy:
   redshift, sizes, data quality) and the DRP table (NASA-Sloan Atlas masses,
   colours, shapes).

4. **MaNGA DR17 velocity maps** for the whole survey, for the hunt for gas
   that rotates against its stars. Each galaxy's MAPS file is a 1-5 MB gzipped
   FITS file with ~60 extensions; we need only nine of them, and they sit in
   the first ~40% of the file. ``manga_io.py`` therefore streams each file,
   stops reading once it has what it needs, and keeps a ~20 kB table per
   galaxy (one row per Voronoi bin). The original file's size is recorded so
   any galaxy can be re-fetched in full. About 8 GB is transferred (~1 hour at
   2 MB/s) and 240 MB kept.

Everything lands in ``data/`` alongside ``manifest.json``, which records each
file's source URL, retrieval time and SHA-256 hash. Files already present are
never re-downloaded; ``--refresh`` forces it. Big downloads resume where they
stopped.

    python fetch_data.py                 # everything (~1.5 h the first time)
    python fetch_data.py --stage things  # one stage: things, refs, manga, survey
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from config import (DATA_DIR, GALAXIES, MANGA_DAP_URL, MANGA_DAPALL, MANGA_DRP_URL,
                    MANGA_DRPALL, MANGA_DIR, MANGA_PAPER, MANIFEST_PATH, REF_DIR,
                    THINGS_DIR, THINGS_PAPER, THINGS_URL, USER_AGENT)

SPARC_URL = "https://astroweb.cwru.edu/SPARC/"
VIZIER_TSV = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source={table}&-out.all&-out.max=unlimited"


@dataclass(frozen=True)
class Source:
    """One file we download, and why."""

    path: Path
    url: str
    what: str
    reference: str
    #: ``True`` if the analysis is built on it; ``False`` if it is a published
    #: result used only to grade our own at the end.
    input: bool = True
    #: Bytes that must appear near the start of the file (cheap check that we
    #: did not save an HTML error page).
    marker: bytes = b""


def _sources() -> dict[str, list[Source]]:
    things = []
    for g in GALAXIES:
        things.append(Source(
            THINGS_DIR / g.cube_file(), THINGS_URL + g.cube_file(),
            f"{g.name}: natural-weighted 21-cm data cube (no primary-beam or residual "
            "flux correction -- the version de Blok et al. 2008 fitted)",
            THINGS_PAPER, marker=b"SIMPLE  ="))
        for m, what in ((0, "integrated HI intensity (moment 0)"),
                        (1, "intensity-weighted mean velocity (moment 1)"),
                        (2, "velocity dispersion (moment 2)")):
            things.append(Source(
                THINGS_DIR / g.moment_file(m), THINGS_URL + g.moment_file(m),
                f"{g.name}: the THINGS team's {what} map, natural weighting",
                THINGS_PAPER, marker=b"SIMPLE  ="))

    refs = [
        Source(REF_DIR / "deblok2008_arxiv0810.2100.tar.gz", "https://arxiv.org/e-print/0810.2100",
               "de Blok et al. 2008 source: tilted-ring table, and the rotation-curve "
               "figures as vector PostScript (fig7 NGC 2403, fig9 NGC 2903, fig12 NGC 3198, "
               "fig15 DDO 154)",
               "de Blok et al. 2008, AJ 136, 2648", input=False, marker=b"\x1f\x8b"),
        Source(REF_DIR / "trachternach2008_arxiv0810.2116.tar.gz", "https://arxiv.org/e-print/0810.2116",
               "Trachternach et al. 2008 source: harmonic decompositions (non-circular "
               "motion amplitudes) and dynamical centres of the THINGS galaxies",
               "Trachternach et al. 2008, AJ 136, 2720", input=False, marker=b"\x1f\x8b"),
        Source(REF_DIR / "oh2011_arxiv1011.0899.tar.gz", "https://arxiv.org/e-print/1011.0899",
               "Oh et al. 2011 source: dark matter density profiles and inner slopes of "
               "the THINGS dwarfs, DDO 154 included",
               "Oh et al. 2011, AJ 141, 193", input=False, marker=b"\x1f\x8b"),
        Source(REF_DIR / "walter2008_table1.tsv", VIZIER_TSV.format(table="J/AJ/136/2563/table1"),
               "THINGS galaxy table: positions and the distances every THINGS paper adopts",
               THINGS_PAPER, marker=b"#"),
        Source(REF_DIR / "walter2008_table4.tsv", VIZIER_TSV.format(table="J/AJ/136/2563/table4"),
               "THINGS global HI spectra from the blanked natural-weighted cubes: checks our "
               "reading of the cubes' velocity axes",
               THINGS_PAPER, input=False, marker=b"#"),
        Source(REF_DIR / "SPARC_Lelli2016c.mrt", SPARC_URL + "SPARC_Lelli2016c.mrt",
               "SPARC galaxy table (distances, inclinations, data references)",
               "Lelli, McGaugh & Schombert 2016, AJ 152, 157", input=False, marker=b"Title:"),
        Source(REF_DIR / "MassModels_Lelli2016c.mrt", SPARC_URL + "MassModels_Lelli2016c.mrt",
               "SPARC rotation curves (project D4's), for the same galaxies from older "
               "radio and optical data", "Lelli, McGaugh & Schombert 2016, AJ 152, 157",
               input=False, marker=b"Title:"),
    ]

    manga = [
        Source(MANGA_DIR / MANGA_DAPALL, MANGA_DAP_URL + MANGA_DAPALL,
               "MaNGA DR17 DAP summary: one row per galaxy and binning scheme -- redshift, "
               "effective radius, data-quality bits", MANGA_PAPER, marker=b"SIMPLE  ="),
        Source(MANGA_DIR / MANGA_DRPALL, MANGA_DRP_URL + MANGA_DRPALL,
               "MaNGA DR17 DRP summary: NASA-Sloan Atlas stellar masses, colours, axis "
               "ratios and position angles", MANGA_PAPER, marker=b"SIMPLE  ="),
    ]
    return {"things": things, "refs": refs, "manga": manga}


# --------------------------------------------------------------------------
# Downloading
# --------------------------------------------------------------------------

_manifest_lock = threading.Lock()


def load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def record(path: Path, url: str, what: str, reference: str, used_as: str, **extra) -> None:
    """Add one file to the provenance manifest (thread-safe)."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            digest.update(chunk)
    entry = {
        "url": url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
        "what": what,
        "reference": reference,
        "used_as": used_as,
        **extra,
    }
    with _manifest_lock:
        manifest = load_manifest()
        manifest["files"][str(path.relative_to(DATA_DIR))] = entry
        tmp = MANIFEST_PATH.with_suffix(".json.part")
        tmp.write_text(json.dumps(manifest, indent=2))
        tmp.replace(MANIFEST_PATH)


def download(url: str, dest: Path, *, marker: bytes = b"", attempts: int = 6,
             quiet: bool = False) -> int:
    """Fetch ``url`` to ``dest``, resuming a partial ``.part`` file; return bytes.

    The file only gets its final name once complete, so an interrupted
    download never looks cached.
    """
    part = dest.with_name(dest.name + ".part")
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        have = part.stat().st_size if part.exists() else 0
        headers = {"User-Agent": USER_AGENT}
        if have:
            headers["Range"] = f"bytes={have}-"
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=120) as response:
                if have and response.status != 206:  # server ignored the range: restart
                    have = 0
                total = have + int(response.headers.get("Content-Length") or 0)
                with open(part, "ab" if have else "wb") as fh:
                    t0, done = time.time(), have
                    while chunk := response.read(1 << 20):
                        fh.write(chunk)
                        done += len(chunk)
                        if not quiet and total > 50e6 and done % (100 << 20) < (1 << 20):
                            rate = (done - have) / max(time.time() - t0, 1e-3) / 1e6
                            print(f"    {dest.name}: {done / 1e6:.0f}/{total / 1e6:.0f} MB "
                                  f"({rate:.1f} MB/s)", flush=True)
            with open(part, "rb") as fh:
                head = fh.read(4096)
            if marker and marker not in head:
                part.unlink()
                raise ValueError(f"{dest.name}: content does not look like the expected file")
            part.replace(dest)
            return dest.stat().st_size
        except (urllib.error.URLError, TimeoutError, ConnectionError, ValueError) as err:
            last_error = err
            if attempt < attempts:
                time.sleep(min(30.0, 3.0 * attempt))
    raise RuntimeError(f"could not download {url}: {last_error}")


def fetch_sources(sources: list[Source], *, refresh: bool = False, workers: int = 1,
                  quiet: bool = False) -> None:
    def one(source: Source) -> None:
        if source.path.exists() and not refresh:
            if not quiet:
                print(f"  cached     {source.path.relative_to(DATA_DIR)}")
            return
        if not quiet:
            print(f"  fetching   {source.path.relative_to(DATA_DIR)}", flush=True)
        size = download(source.url, source.path, marker=source.marker, quiet=quiet)
        record(source.path, source.url, source.what, source.reference,
               "analysis input" if source.input else "published comparison only")
        if not quiet:
            print(f"  done       {source.path.relative_to(DATA_DIR)} ({size / 1e6:.1f} MB)",
                  flush=True)

    if workers == 1:
        for s in sources:
            one(s)
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for future in [pool.submit(one, s) for s in sources]:
                future.result()


def missing(stage: str) -> list[str]:
    return [str(s.path.relative_to(DATA_DIR)) for s in _sources()[stage] if not s.path.exists()]


def verify_cache() -> list[str]:
    """Cached files whose bytes no longer match the recorded hash (big files skipped)."""
    changed = []
    for name, rec in load_manifest().get("files", {}).items():
        path = DATA_DIR / name
        if path.exists() and rec["bytes"] < 200e6:
            if hashlib.sha256(path.read_bytes()).hexdigest() != rec["sha256"]:
                changed.append(name)
    return changed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--stage", action="append",
                        choices=["things", "refs", "manga", "survey"],
                        help="fetch only this stage (repeatable); default: all")
    parser.add_argument("--refresh", action="store_true",
                        help="re-download files even if cached")
    args = parser.parse_args(argv)
    stages = args.stage or ["refs", "things", "manga", "survey"]
    sources = _sources()

    try:
        if "refs" in stages:
            print("Published results (grading only) and SPARC")
            fetch_sources(sources["refs"], refresh=args.refresh)
        if "things" in stages:
            print("THINGS 21-cm cubes and moment maps (~2 GB)")
            fetch_sources(sources["things"], refresh=args.refresh, workers=3)
        if "manga" in stages:
            print("MaNGA DR17 summary catalogues (~220 MB)")
            fetch_sources(sources["manga"], refresh=args.refresh, workers=2)
        if "survey" in stages:
            import manga_io
            print("MaNGA DR17 velocity maps, every galaxy (streamed; ~1 hour)")
            manga_io.fetch_survey(refresh=args.refresh)
    except RuntimeError as err:
        print(f"\n{err}", file=sys.stderr)
        return 1

    changed = verify_cache()
    if changed:
        print(f"\nWARNING: cached files differ from their recorded hashes: {changed}")
    print("\nDone; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
