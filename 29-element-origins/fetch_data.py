"""Download and cache every input of L4 (run once; ~1.1 GB, mostly GALAH).

    python fetch_data.py            # everything not yet cached
    python fetch_data.py --only galah,kilonova

Sources
-------
* GALAH DR4 (Buder et al. 2025), Data Central: the main catalogue, the
  orbit value-added catalogue and the solar zero points.
* GALAH DR4 europium VAC (Kane et al. 2026): DR4's own Eu is unreliable
  (see the README); this catalogue re-derives it.
* Villar et al. (2017): the AT2017gfo photometry table, from the arXiv source.
* SVO Filter Profile Service: the transmission curve of every filter in it.
* Asplund et al. (2009): the solar abundance table, from the arXiv source.
* Bisterzo et al. (2014), VizieR: solar s-process fractions by element.
* JINAbase (Abohalima & Frebel 2018; updated 2021): metal-poor stars.
* APOGEE DR17 (VizieR III/286): [Fe/H], [Mg/Fe] for the cross-check.

Every file is written to a temporary name and renamed when complete, so an
interrupted run never leaves a truncated file that looks cached.
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
import tarfile
import time
import urllib.request

import numpy as np

from config import (APOGEE_DIR, APOGEE_VIZIER, ASPLUND_ARXIV, BISTERZO_VIZIER, EU_DIR,
                    EU_VAC_BASE, EU_VAC_FILES, GALAH_BASE, GALAH_DIR, GALAH_FILES,
                    JINABASE_URL, KN_DIR, LIT_DIR, SVO_URL, VILLAR_ARXIV)

UA = "Mozilla/5.0 (X11; Linux x86_64) astrophysics-atlas-L4"


def _open(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(req, timeout=timeout)


def download(url: str, dest, *, retries: int = 6, quiet: bool = False) -> bool:
    """Fetch ``url`` to ``dest`` (resuming a partial file). Returns True if
    the file was fetched now, False if it was already cached."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return False
    part = dest.with_name(dest.name + ".part")
    for attempt in range(retries):
        try:
            have = part.stat().st_size if part.exists() else 0
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            if have:
                req.add_header("Range", f"bytes={have}-")
            with urllib.request.urlopen(req, timeout=120) as r:
                if have and r.status != 206:      # server ignored the range
                    have = 0
                total = r.headers.get("Content-Length")
                total = int(total) + have if total else None
                with open(part, "ab" if have else "wb") as f:
                    t0, done = time.time(), have
                    while chunk := r.read(1 << 20):
                        f.write(chunk)
                        done += len(chunk)
                        if not quiet and total and total > 50e6:
                            rate = (done - have) / max(time.time() - t0, 1e-3)
                            print(f"\r  {dest.name}: {done / 1e6:7.1f} / {total / 1e6:.1f} MB "
                                  f"({rate / 1e6:.1f} MB/s)", end="", flush=True)
            if total and part.stat().st_size != total:
                raise OSError(f"short read: {part.stat().st_size} of {total} bytes")
            part.rename(dest)
            if not quiet and total and total > 50e6:
                print()
            return True
        except Exception as exc:   # network hiccups: retry with backoff
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} for {dest.name}: {exc}")
            time.sleep(3 * (attempt + 1))
    return False


def _arxiv_source(arxiv_id: str) -> tarfile.TarFile:
    data = _open(f"https://arxiv.org/e-print/{arxiv_id}").read()
    return tarfile.open(fileobj=io.BytesIO(gzip.decompress(data)))


# ---------------------------------------------------------------------------


def fetch_galah():
    for key, name in GALAH_FILES.items():
        if download(f"{GALAH_BASE}/{name}", GALAH_DIR / name):
            print(f"  fetched {name}")
    for name in EU_VAC_FILES:
        if download(f"{EU_VAC_BASE}/{name}", EU_DIR / name):
            print(f"  fetched europium VAC {name}")


def fetch_villar():
    dest = KN_DIR / "villar2017_phot_table.dat"
    if dest.exists():
        return
    tar = _arxiv_source(VILLAR_ARXIV)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(tar.extractfile("phot_table.dat").read())
    print(f"  extracted Villar et al. (2017) Table 3 -> {dest.name}")


def fetch_filters():
    import kn_data

    t = kn_data.read_villar_table()
    ids = sorted(set(t["svo"][t["used"]]))
    n_new = 0
    for svo in ids:
        path = kn_data.filter_path(svo)
        if path.exists():
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        xml = _open(SVO_URL.format(svo)).read().decode()
        rows = re.findall(r"<TR>\s*<TD>([^<]+)</TD>\s*<TD>([^<]+)</TD>\s*</TR>", xml)
        if not rows:
            raise RuntimeError(f"SVO returned no transmission curve for {svo}")
        np.savetxt(path, np.array(rows, float), header=f"{svo} (SVO FPS): wavelength [A], transmission")
        n_new += 1
    if n_new:
        print(f"  fetched {n_new} SVO filter curves ({len(ids)} in use)")


def fetch_asplund():
    """Table 1 of Asplund et al. (2009), parsed from its LaTeX source."""
    dest = LIT_DIR / "asplund2009_table1.csv"
    if dest.exists():
        return
    tex = _arxiv_source(ASPLUND_ARXIV).extractfile("sun.tex").read().decode()
    start = tex.index(r"\label{t:sun}")
    body = tex[start:tex.index(r"\end{tabular}", start)]
    def value(cell):
        # "$8.22 \pm 0.04$", "$[10.93 \pm 0.01]$" (indirect), "-$0.12 \pm 0.04$", "$$"
        m = re.match(r"\s*(-?[\d.]+)", re.sub(r"[$\[\]]", "", cell))
        return f"{float(m.group(1)):.2f}" if m else ""

    out = ["Z,element,photosphere,meteorites"]
    for line in body.splitlines():
        line = line.split("%")[0]
        if not re.match(r"\s*\d+\s*&", line):
            continue
        cells = [c.strip() for c in line.replace("\\\\", "").split("&")]
        for half in (cells[:4], cells[4:8]):
            if len(half) == 4 and half[0].isdigit():
                out.append(f"{int(half[0])},{half[1]},{value(half[2])},{value(half[3])}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(out) + "\n")
    print(f"  parsed Asplund et al. (2009) Table 1: {len(out) - 1} elements")


def fetch_bisterzo():
    dest = LIT_DIR / "bisterzo2014_sprocess.csv"
    if dest.exists():
        return
    from astroquery.vizier import Vizier

    t = Vizier(row_limit=-1, columns=["**"]).get_catalogs(BISTERZO_VIZIER)[0]
    dest.parent.mkdir(parents=True, exist_ok=True)
    t.write(dest, format="ascii.csv", overwrite=True)
    print(f"  fetched Bisterzo et al. (2014) s-process table: {len(t)} rows")


def fetch_jinabase():
    if download(JINABASE_URL, LIT_DIR / "JINAbase_2021.csv"):
        print("  fetched JINAbase (2021 update)")


def fetch_apogee():
    """APOGEE DR17 main-survey stars with S/N > 70: only the columns we use
    (abundances, flags, and Gaia astrometry + APOGEE velocities for a
    kinematic disc selection). VizieR's TAP service is tried first, then its
    ASU interface."""
    dest = APOGEE_DIR / "apogee_dr17_mgfe.fits"
    if dest.exists():
        return
    cols = ['GaiaEDR3', 'RAJ2000', 'DEJ2000', 'Teff', 'logg', '[Fe/H]', 'f_[Fe/H]', '[Mg/Fe]',
            'f_[Mg/Fe]', 'SNR', 'AFlag', 'EFlag', 'Tel', 'plx', 'e_plx', 'pmRA', 'pmDE', 'HRV']
    t = None
    try:
        import pyvo

        svc = pyvo.dal.TAPService("https://tapvizier.cds.unistra.fr/TAPVizieR/tap")
        q = ("SELECT " + ", ".join(f'"{c}"' for c in cols)
             + f' FROM "{APOGEE_VIZIER}" WHERE "SNR" > 70 AND "EFlag" = 0')
        t = svc.run_async(q, maxrec=2_000_000).to_table()
    except Exception as exc:
        print(f"  VizieR TAP failed ({exc}); using the ASU interface")
        from astroquery.vizier import Vizier

        v = Vizier(columns=cols, row_limit=-1, column_filters={"SNR": ">70", "EFlag": "=0"},
                   timeout=1800)
        t = v.get_catalogs(APOGEE_VIZIER)[0]
    dest.parent.mkdir(parents=True, exist_ok=True)
    for c in t.colnames:
        t[c].description = None
    t.write(dest, overwrite=True)
    print(f"  fetched APOGEE DR17: {len(t):,} stars (S/N > 70, main survey)")


STEPS = {
    "galah": fetch_galah,
    "kilonova": lambda: (fetch_villar(), fetch_filters()),
    "literature": lambda: (fetch_asplund(), fetch_bisterzo(), fetch_jinabase()),
    "apogee": fetch_apogee,
}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--only", help="comma-separated subset of: " + ", ".join(STEPS))
    args = p.parse_args(argv)
    steps = args.only.split(",") if args.only else list(STEPS)
    for s in steps:
        print(f"[{s}]")
        STEPS[s]()
    print("all inputs cached.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
