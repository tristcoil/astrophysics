"""Download and cache every input (about 3.5 MB in all).

    python fetch_data.py          # fetch what is missing
    python fetch_data.py --force  # fetch everything again

What comes from where:

* Saikawa & Shirai (2018): g*rho(T), g*s(T) and w(T) of the Standard Model, from the authors' page.
* Laine, Meyer & Schroeder (2015, updated 2019): the Standard Model equation of state, from the
  authors' page.
* Borsanyi et al. (2014): the 2+1 flavour lattice QCD equation of state (arXiv ancillary file).
* Escudero (2020): neutrino decoupling, T_gamma/T_nu, g* and g_s from 10 MeV to 1 keV (arXiv
  ancillary file).
* HotQCD (Bazavov et al. 2014): the 2+1 flavour lattice equation of state, Table I, read from
  the paper's arXiv source.
* Borsanyi et al. (2016): the 2+1+1 flavour trace-anomaly parametrisation and their spline table
  of g_rho(T) and g_rho/g_s for the whole Standard Model, read from the paper's arXiv source.

The PDG particle list needs no download: it ships with the `particle` package.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import config as C  # noqa: E402

UA = "Mozilla/5.0 (X11; Linux x86_64) astrophysics-atlas-P2"

#: (local name, url, what it is)
FILES = [
    ("ss2018_standardmodel.dat", "https://member.ipmu.jp/satoshi.shirai/standardmodel2018.dat",
     "Saikawa & Shirai 2018: T, g*rho, err, g*s, err"),
    ("ss2018_eos_w.dat", "https://member.ipmu.jp/satoshi.shirai/eos_w.dat",
     "Saikawa & Shirai 2018 (added 2024-07-19): T, w, err"),
    ("lms2015_standardmodel.dat", "http://www.laine.itp.unibe.ch/eos15/standardmodel.dat",
     "Laine, Meyer & Schroeder 2015 (Dec 2019 update): T/MeV, p/T^4, e/T^4, s/T^3, c/T^3, w, cs2, g_eff, h_eff, i_eff"),
    ("wb2014_eos.dat", "https://arxiv.org/src/1309.5258/anc/WB-EoS.dat",
     "Borsanyi et al. 2014 (PLB 730, 99): 2+1 flavour lattice EoS table"),
    ("escudero2020_sm_thermo.dat", "https://arxiv.org/src/2001.04466/anc/SM_Thermo.dat",
     "Escudero 2020 (JCAP 05, 048): neutrino decoupling in the SM"),
    ("ejls2025_rate_coefficients.dat",
     "https://raw.githubusercontent.com/MiguelEA/nudec_BSM/9f1824ec34031642feb0f626a162931bc48bb648/"
     "BasicModules_data/rate_coefficients_neutrinos.dat",
     "Escudero, Jackson, Laine & Sandner 2025 (arXiv:2511.04747): electron-mass dependent "
     "neutrino energy-transfer coefficients vs m_e/T (nudec_BSM v2, pinned commit)"),
]

#: (arXiv id, members to keep, what it is). Only the named members are extracted.
ARXIV_SOURCES = [
    ("1407.6387", ["tab_fin.tex"], "HotQCD 2014 (PRD 90, 094503): Table I, continuum EoS"),
    ("1606.07494", ["som_eos.tex"], "Borsanyi et al. 2016 (Nature 539, 69): EoS supplement"),
]


def _get(url: str, *, tries: int = 4, timeout: float = 90.0) -> bytes:
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # network hiccups: retry with a pause
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError(f"could not fetch {url}: {last}")


def _sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _load_manifest() -> dict:
    if C.MANIFEST_PATH.exists():
        return json.loads(C.MANIFEST_PATH.read_text())
    return {}


def fetch_all(force: bool = False) -> dict:
    C.RAW_DIR.mkdir(parents=True, exist_ok=True)
    C.LIT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = _load_manifest()
    today = dt.date.today().isoformat()

    for name, url, what in FILES:
        path = C.RAW_DIR / name
        if path.exists() and not force:
            continue
        print(f"  fetching {url}")
        b = _get(url)
        if len(b) < 1000 or b.lstrip()[:1] == b"<":
            raise RuntimeError(f"{url} did not return a data table ({len(b)} bytes)")
        path.write_bytes(b)
        manifest[f"raw/{name}"] = {"url": url, "what": what, "bytes": len(b),
                                   "sha256": _sha256(b), "fetched": today}

    for arxiv_id, members, what in ARXIV_SOURCES:
        outdir = C.LIT_DIR / arxiv_id
        if all((outdir / m).exists() for m in members) and not force:
            continue
        url = f"https://arxiv.org/e-print/{arxiv_id}"
        print(f"  fetching {url}")
        b = _get(url)
        outdir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(fileobj=io.BytesIO(b), mode="r:gz") as tf:
            for m in members:
                info = tf.getmember(m)  # exact name only: nothing else is unpacked
                if not info.isfile():
                    raise RuntimeError(f"{arxiv_id}: {m} is not a regular file")
                data = tf.extractfile(info).read()
                (outdir / m).write_bytes(data)
                manifest[f"literature/{arxiv_id}/{m}"] = {
                    "url": url, "what": what, "bytes": len(data), "sha256": _sha256(data),
                    "source_tarball_bytes": len(b), "fetched": today}

    try:
        import particle

        manifest["particle_package"] = {"version": particle.__version__,
                                        "table": list(particle.Particle.table_names()),
                                        "what": "PDG particle data bundled with scikit-hep/particle"}
    except ImportError:
        print("  note: `pip install particle` (the PDG table) is needed for the analysis")

    C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    m = fetch_all(force=args.force)
    total = sum(v.get("source_tarball_bytes", v.get("bytes", 0)) for v in m.values() if isinstance(v, dict))
    print(f"  {len(m)} entries in the manifest, {total / 1e6:.1f} MB downloaded in all")
    return 0


if __name__ == "__main__":
    sys.exit(main())
