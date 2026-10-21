"""Download N6's one new data set and fingerprint every input from earlier projects.

    python fetch_data.py

* GW250114 posterior samples (LIGO-Virgo-KAGRA, NRSur7dq4 waveform, the
  discovery paper's preferred analysis; Zenodo record 16877102, 27 MB) and
  GWOSC's parameter table for the event (a few kB).
* SHA-256 fingerprints of the files N6 reads from projects H1, I1-I5, C1 and
  D7, so a later change in their results is noticed (``run.py`` warns).

Everything is recorded in ``data/manifest.json`` (URL, retrieval time, size,
SHA-256). Re-running skips files that are already present and intact.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
import time
import urllib.request

import config as C


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest, *, marker: bytes = b"", tries: int = 4) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": C.USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            head = tmp.read_bytes()[:64]
            if marker and marker not in head:
                raise ValueError(f"unexpected content at {url}: {head[:40]!r}")
            tmp.replace(dest)
            return
        except Exception as exc:  # noqa: BLE001 -- network errors of every kind: retry
            print(f"  attempt {attempt} failed: {exc}")
            time.sleep(5 * attempt)
    raise RuntimeError(f"could not download {url}")


def main() -> dict:
    manifest = json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {}
    downloads = manifest.setdefault("downloads", {})
    jobs = [
        (C.GW250114["samples_file"], C.GW250114["samples_url"], b"\x89HDF",
         "GW250114 posterior samples, NRSur7dq4 (the discovery paper's preferred waveform). "
         + C.GW250114["paper"] + "; data: " + C.GW250114["zenodo_record"]),
        ("gw250114_gwosc_parameters.json", C.GW250114["gwosc_parameters_url"], b"{",
         "GWOSC summary parameters of GW250114 for every waveform model (v1, discovery papers)"),
    ]
    for name, url, marker, what in jobs:
        dest = C.RAW_DIR / name
        rec = downloads.get(name)
        if dest.exists() and rec and rec.get("sha256") == sha256(dest):
            print(f"have {name}")
            continue
        print(f"downloading {name} ...")
        t0 = time.time()
        download(url, dest, marker=marker)
        downloads[name] = {"url": url, "what": what, "bytes": dest.stat().st_size,
                           "sha256": sha256(dest), "retrieved_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
                           "seconds": round(time.time() - t0, 1)}
        print(f"  {dest.stat().st_size / 1e6:.1f} MB in {time.time() - t0:.0f} s")

    inputs = {}
    for label, path in C.INPUT_FILES.items():
        if path.exists():
            inputs[label] = {"path": str(path.relative_to(C.ATLAS_DIR)), "bytes": path.stat().st_size,
                             "sha256": sha256(path)}
        else:
            inputs[label] = {"path": str(path.relative_to(C.ATLAS_DIR)), "missing": True}
            print(f"note: {label} not found ({path}); the analysis will use fallback values")
    manifest["inputs_from_earlier_projects"] = inputs
    manifest["fingerprinted_utc"] = dt.datetime.now(dt.UTC).isoformat(timespec="seconds")
    C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1))
    print(f"manifest -> {C.MANIFEST_PATH}")
    return manifest


def check_inputs(verbose: bool = True) -> list:
    """Compare the earlier projects' files with the fingerprints taken at fetch time."""
    if not C.MANIFEST_PATH.exists():
        return ["no manifest: run fetch_data.py"]
    rec = json.loads(C.MANIFEST_PATH.read_text()).get("inputs_from_earlier_projects", {})
    changed = []
    for label, path in C.INPUT_FILES.items():
        r = rec.get(label)
        if r is None or r.get("missing"):
            continue
        if not path.exists():
            changed.append(f"{label}: missing now")
        elif sha256(path) != r["sha256"]:
            changed.append(f"{label}: changed since fetch_data.py fingerprinted it")
    if verbose:
        for c_ in changed:
            print("WARNING:", c_)
    return changed


if __name__ == "__main__":
    main()
    sys.exit(0)
