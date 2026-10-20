"""MaNGA DR17: the survey catalogues, and every galaxy's velocity maps, streamed.

Each MaNGA galaxy comes with a DAP "MAPS" file: a gzipped FITS file of ~60
extensions (fluxes, velocities and dispersions of 35 emission lines, 46
spectral indices, ...), 1-5 MB compressed. For stars-vs-gas rotation we need
nine of them: where each Voronoi bin is, its S/N, the stellar velocity and
dispersion, and the H-alpha velocity and flux, with their errors and masks.
Those sit in the first ~40% of the file, so we decompress the download as it
arrives, parse FITS headers on the fly, and hang up once the last needed
extension has been read. One small table per galaxy is kept (one row per bin,
~20 kB), so the survey costs ~8 GB of transfer (35-40% of the files) and
~240 MB of disk.

VOR10 binning: spaxels are grouped (Voronoi tessellation) until each bin
reaches S/N 10 in the continuum; stars and gas are then measured on the same
bins. Every spaxel of a bin carries the bin's values, so one spaxel per bin
is enough.
"""

from __future__ import annotations

import hashlib
import json
import threading
import time
import urllib.error
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits

from config import (MANGA_BINS_DIR, MANGA_DAP_URL, MANGA_DAPALL, MANGA_DAPTYPE, MANGA_DIR,
                    MANGA_DRPALL, MANGA_PAPER, MANGA_WORKERS, USER_AGENT)

#: Extensions we keep, in file order. The last one ends the download.
NEEDED = ("BINID", "BIN_LWSKYCOO", "BIN_LWELLCOO", "BIN_SNR",
          "STELLAR_VEL", "STELLAR_VEL_IVAR", "STELLAR_VEL_MASK",
          "STELLAR_SIGMA", "STELLAR_SIGMA_IVAR", "STELLAR_SIGMA_MASK", "STELLAR_SIGMACORR",
          "EMLINE_GFLUX", "EMLINE_GFLUX_IVAR", "EMLINE_GFLUX_MASK",
          "EMLINE_GVEL", "EMLINE_GVEL_IVAR", "EMLINE_GVEL_MASK")
HALPHA = "Ha-6564"

SURVEY_INDEX = MANGA_DIR / "survey_index.jsonl"
SURVEY_FAILURES = MANGA_DIR / "survey_failures.jsonl"

_DTYPES = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}


def maps_url(plate: int, ifu: int) -> str:
    return (f"{MANGA_DAP_URL}{MANGA_DAPTYPE}/{plate}/{ifu}/"
            f"manga-{plate}-{ifu}-MAPS-{MANGA_DAPTYPE}.fits.gz")


# --------------------------------------------------------------------------
# Streaming FITS reader
# --------------------------------------------------------------------------

class _Stream:
    """Decompressed bytes of a gzipped HTTP response, read on demand."""

    def __init__(self, response):
        self.response = response
        self.inflate = zlib.decompressobj(16 + zlib.MAX_WBITS)
        self.buffer = bytearray()
        self.compressed_read = 0
        self.digest = hashlib.sha256()

    def take(self, n: int) -> bytes:
        while len(self.buffer) < n:
            chunk = self.response.read(1 << 16)
            if not chunk:
                self.buffer += self.inflate.flush()
                if len(self.buffer) < n:
                    raise EOFError("stream ended inside a FITS unit")
                break
            self.compressed_read += len(chunk)
            self.digest.update(chunk)
            self.buffer += self.inflate.decompress(chunk)
        out = bytes(self.buffer[:n])
        del self.buffer[:n]
        return out


def _read_hdus(stream: _Stream, wanted: tuple[str, ...]) -> dict[str, tuple[fits.Header, np.ndarray]]:
    """Read FITS units in order until every ``wanted`` extension has been seen."""
    found: dict[str, tuple[fits.Header, np.ndarray]] = {}
    first = True
    while len(found) < len(wanted):
        raw = b""
        while True:
            block = stream.take(2880)
            raw += block
            if any(block[i:i + 8] == b"END     " for i in range(0, 2880, 80)):
                break
        header = fits.Header.fromstring(raw.decode("ascii"))
        naxis = header.get("NAXIS", 0)
        shape = [header[f"NAXIS{k}"] for k in range(naxis, 0, -1)]
        nbytes = abs(header["BITPIX"]) // 8 * int(np.prod(shape)) if naxis else 0
        nbytes += abs(header["BITPIX"]) // 8 * header.get("PCOUNT", 0) * header.get("GCOUNT", 1) \
            if not first else 0
        padded = (nbytes + 2879) // 2880 * 2880
        payload = stream.take(padded) if padded else b""
        name = "PRIMARY" if first else header.get("EXTNAME", "").strip()
        first = False
        if name in wanted:
            data = np.frombuffer(payload[:nbytes], dtype=_DTYPES[header["BITPIX"]]).reshape(shape)
            found[name] = (header, data.astype(data.dtype.newbyteorder("=")))
    return found


def _channel(header: fits.Header, name: str) -> int:
    """Index of the channel called ``name`` in a multi-channel DAP extension."""
    for key, value in header.items():
        if key.startswith("C") and key[1:].isdigit() and str(value).strip() == name:
            return int(key[1:]) - 1
    raise KeyError(f"channel {name!r} not found")


def bins_from_maps(hdus: dict[str, tuple[fits.Header, np.ndarray]]) -> dict[str, np.ndarray]:
    """One row per Voronoi bin: position, S/N, stellar and H-alpha kinematics."""
    binid = hdus["BINID"][1]
    # Channel 1 = stellar-continuum bins, channel 3 = emission-line bins. For VOR10
    # they coincide; check, and use the stellar ones.
    bins_star, bins_gas = binid[1], binid[3]
    if not np.array_equal(bins_star, bins_gas):
        raise ValueError("stellar and emission-line bins differ")
    ids, first = np.unique(bins_star.ravel(), return_index=True)
    keep = ids >= 0
    ids, first = ids[keep], first[keep]
    npix = np.bincount(bins_star.ravel()[bins_star.ravel() >= 0])[ids]

    def per_bin(name: str, channel: int | None = None) -> np.ndarray:
        data = hdus[name][1]
        if channel is not None:
            data = data[channel]
        return data.ravel()[first]

    sky_h, ell_h = hdus["BIN_LWSKYCOO"][0], hdus["BIN_LWELLCOO"][0]
    ha = _channel(hdus["EMLINE_GVEL"][0], HALPHA)
    if _channel(hdus["EMLINE_GFLUX"][0], HALPHA) != ha:
        raise ValueError("H-alpha channel differs between flux and velocity extensions")
    return {
        "bin_id": ids.astype(np.int32),
        "n_spaxels": npix.astype(np.int16),
        "x_arcsec": per_bin("BIN_LWSKYCOO", _channel(sky_h, "Lum. weighted on-sky X")).astype(np.float32),
        "y_arcsec": per_bin("BIN_LWSKYCOO", _channel(sky_h, "Lum. weighted on-sky Y")).astype(np.float32),
        "r_re": per_bin("BIN_LWELLCOO", _channel(ell_h, "R/Re")).astype(np.float32),
        "snr": per_bin("BIN_SNR").astype(np.float32),
        "star_vel": per_bin("STELLAR_VEL").astype(np.float32),
        "star_vel_ivar": per_bin("STELLAR_VEL_IVAR").astype(np.float32),
        "star_vel_mask": per_bin("STELLAR_VEL_MASK").astype(np.int32),
        "star_sigma": per_bin("STELLAR_SIGMA").astype(np.float32),
        "star_sigma_ivar": per_bin("STELLAR_SIGMA_IVAR").astype(np.float32),
        "star_sigma_mask": per_bin("STELLAR_SIGMA_MASK").astype(np.int32),
        # Quadrature corrections for instrumental resolution: channel 0 from the
        # template/galaxy resolution difference, channel 1 from a fit.
        "star_sigmacorr0": per_bin("STELLAR_SIGMACORR", 0).astype(np.float32),
        "star_sigmacorr1": per_bin("STELLAR_SIGMACORR", 1).astype(np.float32),
        "ha_vel": per_bin("EMLINE_GVEL", ha).astype(np.float32),
        "ha_vel_ivar": per_bin("EMLINE_GVEL_IVAR", ha).astype(np.float32),
        "ha_vel_mask": per_bin("EMLINE_GVEL_MASK", ha).astype(np.int32),
        "ha_flux": per_bin("EMLINE_GFLUX", ha).astype(np.float32),
        "ha_flux_ivar": per_bin("EMLINE_GFLUX_IVAR", ha).astype(np.float32),
        "ha_flux_mask": per_bin("EMLINE_GFLUX_MASK", ha).astype(np.int32),
    }


def fetch_galaxy_bins(plate: int, ifu: int, *, attempts: int = 4) -> dict:
    """Stream one MAPS file, keep its per-bin table, return provenance."""
    url = maps_url(plate, ifu)
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=120) as response:
                size = int(response.headers.get("Content-Length") or 0)
                stream = _Stream(response)
                hdus = _read_hdus(stream, NEEDED)
            table = bins_from_maps(hdus)
            out = MANGA_BINS_DIR / f"{plate}-{ifu}.npz"
            tmp = out.with_name(out.name + ".part.npz")
            np.savez_compressed(tmp, **table)
            tmp.replace(out)
            return {"plateifu": f"{plate}-{ifu}", "url": url, "file_bytes": size,
                    "bytes_read": stream.compressed_read,
                    "sha256_of_bytes_read": stream.digest.hexdigest(),
                    "n_bins": int(len(table["bin_id"])),
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        except (urllib.error.URLError, TimeoutError, ConnectionError, EOFError, zlib.error) as err:
            last = err
            if isinstance(err, urllib.error.HTTPError) and err.code == 404:
                break
            time.sleep(3.0 * attempt)
    raise RuntimeError(f"{plate}-{ifu}: {last}")


def load_bins(plateifu: str) -> dict[str, np.ndarray]:
    with np.load(MANGA_BINS_DIR / f"{plateifu}.npz") as f:
        return {k: f[k] for k in f.files}


# --------------------------------------------------------------------------
# Catalogues and the survey sample
# --------------------------------------------------------------------------

#: Critical-failure bits of DRP3QUAL (bit 30) and DAPQUAL (bit 30).
CRITICAL_BIT = 1 << 30


def read_dapall() -> np.ndarray:
    with fits.open(MANGA_DIR / MANGA_DAPALL, memmap=True) as hdul:
        ext = MANGA_DAPTYPE if MANGA_DAPTYPE in hdul else 1
        return np.array(hdul[ext].data)


def read_drpall() -> np.ndarray:
    with fits.open(MANGA_DIR / MANGA_DRPALL, memmap=True) as hdul:
        return np.array(hdul[1].data)


def survey_targets() -> list[tuple[int, int]]:
    """Every galaxy the DAP finished without a critical failure, once each.

    Repeat observations of the same galaxy (same MANGAID) keep the one with
    the highest median S/N.
    """
    dap = read_dapall()
    ok = dap["DAPDONE"].astype(bool)
    ok &= (dap["DAPQUAL"] & CRITICAL_BIT) == 0
    ok &= (dap["DRP3QUAL"] & CRITICAL_BIT) == 0
    if "DAPTYPE" in dap.dtype.names:
        ok &= np.char.strip(dap["DAPTYPE"].astype(str)) == MANGA_DAPTYPE
    rows = dap[ok]
    best: dict[str, tuple[float, int, int]] = {}
    for row in rows:
        mid = row["MANGAID"].decode() if isinstance(row["MANGAID"], bytes) else str(row["MANGAID"])
        mid = mid.strip()
        snr = float(row["SNR_MED"][1]) if np.ndim(row["SNR_MED"]) else float(row["SNR_MED"])
        plate, ifu = int(row["PLATE"]), int(row["IFUDESIGN"])
        if mid not in best or snr > best[mid][0]:
            best[mid] = (snr, plate, ifu)
    return sorted((p, i) for _, p, i in best.values())


def fetch_survey(*, refresh: bool = False, workers: int = MANGA_WORKERS, limit: int | None = None) -> None:
    """Stream every target's MAPS file; resumable (finished galaxies are skipped)."""
    import fetch_data  # for the manifest helpers

    targets = survey_targets()
    if limit:
        targets = targets[:limit]
    done = {p.stem for p in MANGA_BINS_DIR.glob("*.npz") if not p.name.endswith(".part.npz")}
    todo = [(p, i) for p, i in targets if refresh or f"{p}-{i}" not in done]
    print(f"  {len(targets):,} galaxies in DR17 with a finished DAP run; "
          f"{len(targets) - len(todo):,} cached, {len(todo):,} to stream", flush=True)
    if not todo:
        return

    lock = threading.Lock()
    t0, n_done, n_bytes, n_fail = time.time(), 0, 0, 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch_galaxy_bins, p, i): (p, i) for p, i in todo}
        for fut in as_completed(futures):
            p, i = futures[fut]
            try:
                rec = fut.result()
                with lock, open(SURVEY_INDEX, "a") as fh:
                    fh.write(json.dumps(rec) + "\n")
                n_done += 1
                n_bytes += rec["bytes_read"]
            except Exception as err:  # keep going; failures are listed for a re-run
                n_fail += 1
                with lock, open(SURVEY_FAILURES, "a") as fh:
                    fh.write(json.dumps({"plateifu": f"{p}-{i}", "error": str(err)}) + "\n")
            if (n_done + n_fail) % 250 == 0:
                dt = time.time() - t0
                print(f"    {n_done + n_fail:,}/{len(todo):,} galaxies, {n_bytes / 1e9:.2f} GB read, "
                      f"{n_bytes / 1e6 / dt:.1f} MB/s, {n_fail} failed", flush=True)

    n_cached = len(list(MANGA_BINS_DIR.glob("*-*.npz")))
    fetch_data.record(SURVEY_INDEX, MANGA_DAP_URL + MANGA_DAPTYPE + "/",
                      f"Index of the streamed MaNGA DR17 {MANGA_DAPTYPE} MAPS files: one line per "
                      f"galaxy with its URL, full size, bytes read and their SHA-256; "
                      f"{n_cached:,} per-bin tables in manga/bins/",
                      MANGA_PAPER, "analysis input")
    print(f"  streamed {n_done:,} galaxies ({n_bytes / 1e9:.2f} GB) in {(time.time() - t0) / 60:.0f} min; "
          f"{n_fail} failed (see {SURVEY_FAILURES.name})")
