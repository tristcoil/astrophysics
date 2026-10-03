"""Download and cache everything the analysis needs.

Five things are fetched, all free and public, all from the LIGO/Virgo/KAGRA
open-data services:

1. **Strain data** for the two LIGO detectors, Hanford (H1) and Livingston (L1):
   4096 seconds sampled 4096 times per second, centred on the event, from the
   Gravitational Wave Open Science Center (GWOSC). *Strain* is the detector's
   calibrated output: the fractional change in the length difference of its two
   4 km arms. This is the only large download, ~130 MB per detector.

2. **Event parameters** from the GWOSC Event API: the GWTC-1 catalogue's
   published masses, distance and significance, plus the exact file URLs. The
   URLs are taken from this record rather than hard-coded.

3. **The data behind Figure 1 of the discovery paper**, which GWOSC hosts. This
   is the collaboration's own filtered strain, and we check our processing of the
   raw data against it.

4. **A waveform template** for GW150914 that GWOSC provides for its tutorials:
   an SEOBNRv2 waveform, the model family the discovery search itself used.

5. **Posterior samples** from the GWTC-1 parameter estimation: 8350 draws of the
   masses and sky position consistent with the signal. They give the published
   detector-frame chirp mass to compare our own estimate against.

Everything is cached under ``data/``. ``data/manifest.json`` records where each
file came from, when, its size and its SHA-256, so a result can always be traced
to the exact bytes it was computed from. Nothing already present is re-fetched.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from config import (
    DATA_DIR,
    EVENT,
    EVENT_JSON_PATH,
    MANIFEST_PATH,
    POSTERIOR_PATH,
    PRL_DIR,
    STRAIN_DIR,
    TEMPLATE_PATH,
)

EVENT_API_URL = (
    f"https://gwosc.org/eventapi/json/{EVENT.catalog}/{EVENT.name}/v{EVENT.version}/"
)
PRL_FIGURE1_URL = "https://gwosc.org/s/events/GW150914/P150914/"
PRL_FIGURE1_FILES = [
    f"fig1-{kind}-{ifo}.txt"
    for kind in ("observed", "waveform", "residual")
    for ifo in ("H", "L")
]
TEMPLATE_URL = "https://gwosc.org/s/events/GW150914/GW150914_4_template.hdf5"

USER_AGENT = "astrophysics-project-atlas/0.1 (hobbyist; caches every download)"


# --------------------------------------------------------------------------
# Downloading
# --------------------------------------------------------------------------

def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _save_manifest(manifest: dict) -> None:
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, *, quiet: bool = False, refresh: bool = False) -> Path:
    """Fetch ``url`` to ``dest`` unless an intact copy is already cached.

    A cached file counts as intact if its size matches the manifest. Downloads go
    to a temporary name first and are renamed only once complete, so an
    interrupted download can never be mistaken for a finished one.
    """
    manifest = _load_manifest()
    key = str(dest.relative_to(DATA_DIR))
    entry = manifest["files"].get(key)
    if dest.exists() and not refresh and entry and entry.get("bytes") == dest.stat().st_size:
        if not quiet:
            print(f"  cached   {key}")
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        expected = int(response.headers.get("Content-Length") or 0)
        if not quiet:
            size = f"{expected / 1e6:.1f} MB" if expected else "unknown size"
            print(f"  fetching {key} ({size})")
        done, next_report = 0, 0.1
        with open(partial, "wb") as out:
            while chunk := response.read(1 << 20):
                out.write(chunk)
                done += len(chunk)
                if not quiet and expected > 20e6 and done / expected >= next_report:
                    print(f"           {done / expected:4.0%}")
                    next_report += 0.1
    if expected and partial.stat().st_size != expected:
        partial.unlink()
        raise RuntimeError(f"truncated download of {url}: "
                           f"{partial.stat().st_size} of {expected} bytes")
    partial.replace(dest)

    manifest["files"][key] = {
        "url": url,
        "bytes": dest.stat().st_size,
        "sha256": _sha256(dest),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    _save_manifest(manifest)
    return dest


def fetch_event(*, quiet: bool = False, refresh: bool = False) -> dict:
    """The GWOSC Event API record for the event: parameters and file URLs."""
    download(EVENT_API_URL, EVENT_JSON_PATH, quiet=quiet, refresh=refresh)
    return load_event()


def strain_url(event_record: dict, detector: str) -> str:
    """Pick the HDF5 strain file matching the configured rate and duration."""
    for s in event_record["strain"]:
        if (s["detector"] == detector and s["format"] == "hdf5"
                and s["sampling_rate"] == EVENT.sample_rate
                and s["duration"] == EVENT.duration):
            return s["url"]
    raise KeyError(f"no {EVENT.duration} s, {EVENT.sample_rate} Hz HDF5 file for {detector}")


def strain_path(detector: str) -> Path:
    return STRAIN_DIR / f"{detector}_{EVENT.sample_rate}Hz_{EVENT.duration}s.hdf5"


def fetch_all(*, quiet: bool = False, refresh: bool = False) -> None:
    if not quiet:
        print(f"GWOSC event record for {EVENT.name}")
    record = fetch_event(quiet=quiet, refresh=refresh)

    if not quiet:
        print(f"\nStrain: {EVENT.duration} s at {EVENT.sample_rate} Hz per detector")
    for det in EVENT.detectors:
        download(strain_url(record, det), strain_path(det), quiet=quiet, refresh=refresh)

    if not quiet:
        print("\nData behind Figure 1 of the discovery paper")
    for name in PRL_FIGURE1_FILES:
        download(PRL_FIGURE1_URL + name, PRL_DIR / name, quiet=quiet, refresh=refresh)

    if not quiet:
        print("\nWaveform template and GWTC-1 posterior samples")
    download(TEMPLATE_URL, TEMPLATE_PATH, quiet=quiet, refresh=refresh)
    download(posterior_url(record), POSTERIOR_PATH, quiet=quiet, refresh=refresh)


def posterior_url(event_record: dict) -> str:
    """URL of the preferred parameter-estimation release, from the event record."""
    for p in event_record["parameters"].values():
        if p.get("is_preferred") and p.get("pipeline_type") == "pe" and p.get("data_url"):
            return p["data_url"]
    raise KeyError("event record lists no preferred posterior-sample release")


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_event() -> dict:
    """The event's catalogue record (one event, unwrapped from the API envelope)."""
    if not EVENT_JSON_PATH.exists():
        raise FileNotFoundError(f"{EVENT_JSON_PATH} not found -- run `python fetch_data.py` first")
    events = json.loads(EVENT_JSON_PATH.read_text())["events"]
    (record,) = events.values()
    return record


@dataclass
class Strain:
    """One detector's strain time series with its data-quality flags.

    ``values`` is dimensionless strain h(t), sampled uniformly from ``t0`` (GPS
    seconds) every ``dt`` seconds. The quality masks are one integer per second
    whose bits say whether that second passed each data-quality category and
    whether any hardware injection -- a fake signal deliberately driven into the
    mirrors to test the pipelines -- was active.
    """

    detector: str
    t0: float
    dt: float
    values: np.ndarray
    dq_mask: np.ndarray
    dq_names: list[str]
    inj_mask: np.ndarray
    inj_names: list[str]
    meta: dict = field(default_factory=dict)

    @property
    def sample_rate(self) -> float:
        return 1.0 / self.dt

    @property
    def times(self) -> np.ndarray:
        return self.t0 + self.dt * np.arange(len(self.values))

    def flag(self, name: str, *, injections: bool = False) -> np.ndarray:
        """Per-second boolean for one named quality or injection bit."""
        names = self.inj_names if injections else self.dq_names
        mask = self.inj_mask if injections else self.dq_mask
        return (mask >> names.index(name)) & 1 == 1


def load_strain(detector: str) -> Strain:
    """Read a GWOSC HDF5 strain file.

    GWOSC's format keeps the strain in ``strain/Strain``, with its start time and
    sample spacing as attributes, and the 1 Hz quality bitmasks under
    ``quality/``. Each bit's name is stored alongside the mask, so bits are
    looked up by name rather than by hard-coded position.
    """
    import h5py

    path = strain_path(detector)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")

    def names(ds) -> list[str]:
        return [n.decode() if isinstance(n, bytes) else str(n) for n in ds[()]]

    with h5py.File(path, "r") as f:
        ds = f["strain/Strain"]
        values = ds[()].astype(np.float64)
        t0 = float(ds.attrs["Xstart"])
        dt = float(ds.attrs["Xspacing"])
        meta = {k: (v.decode() if isinstance(v, bytes) else v.item() if hasattr(v, "item") else v)
                for k, v in ((k, f["meta"][k][()]) for k in f["meta"])}
        strain = Strain(
            detector=detector, t0=t0, dt=dt, values=values,
            dq_mask=f["quality/simple/DQmask"][()].astype(np.int64),
            dq_names=names(f["quality/simple/DQShortnames"]),
            inj_mask=f["quality/injections/Injmask"][()].astype(np.int64),
            inj_names=names(f["quality/injections/InjShortnames"]),
            meta=meta,
        )
    if not np.isclose(strain.sample_rate, EVENT.sample_rate):
        raise ValueError(f"{path} is {strain.sample_rate} Hz, expected {EVENT.sample_rate}")
    return strain


def load_template() -> dict:
    """GWOSC's GW150914 template: plus and cross polarisations plus metadata.

    The two rows are the waveform's two polarisations, h+ and hx, which are the
    same chirp 90 degrees out of phase. Matched filtering uses both so that the
    unknown phase of the real signal can be maximised over analytically.
    """
    import h5py

    with h5py.File(TEMPLATE_PATH, "r") as f:
        plus, cross = f["template"][()].astype(np.float64)
        meta = {k: (v.item() if hasattr(v, "item") else v) for k, v in f["meta"].attrs.items()}
    return {"plus": plus, "cross": cross, "meta": meta, "dt": 1.0 / float(meta["fs"])}


def load_posterior(key: str = "Overall_posterior") -> np.ndarray:
    """GWTC-1 posterior samples (a structured array, one row per draw)."""
    import h5py

    with h5py.File(POSTERIOR_PATH, "r") as f:
        return f[key][()]


def load_prl_figure1() -> dict:
    """The discovery paper's Figure 1 time series, keyed like ``observed-H``.

    Each file holds time in seconds after 2015-09-14 09:50:45 UTC and strain in
    units of 1e-21, band-passed 35-350 Hz with the instrumental lines notched,
    exactly as plotted in the paper. They are sampled at 16384 Hz.
    """
    out = {}
    for name in PRL_FIGURE1_FILES:
        arr = np.loadtxt(PRL_DIR / name)
        out[name.removeprefix("fig1-").removesuffix(".txt")] = {
            "t": arr[:, 0], "strain": arr[:, 1] * 1e-21,
        }
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download even files that are already cached")
    args = parser.parse_args(argv)

    fetch_all(refresh=args.refresh)
    total = sum(e["bytes"] for e in _load_manifest()["files"].values())
    print(f"\nDone: {total / 1e6:.0f} MB cached in data/. Analysis can now run offline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
