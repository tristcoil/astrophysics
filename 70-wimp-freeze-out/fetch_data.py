"""Download and cache every input (about 25 MB in all).

    python fetch_data.py          # fetch what is missing
    python fetch_data.py --force  # fetch everything again

What comes from where (all public, no accounts):

* **The Standard Model's equation of state**: g*rho(T) and g*s(T) from Saikawa & Shirai (2018,
  2020) and Laine & Meyer (2015, Dec 2019 update), from the authors' pages, and
  Escudero's (2020) neutrino-decoupling history (arXiv ancillary file). P2 (68-thermal-history)
  uses the same tables; when its cache is present the files are hard-linked, not fetched.
* **Answer keys for the relic cross section**: Saikawa & Shirai (2020) for four particle types
  and s/p-wave, and Steigman, Dasgupta & Beacom (2012)'s ancillary curve.
* **Gamma rays from dwarf galaxies**: McDaniel et al. (2024)'s data release on Figshare (the
  per-dwarf, per-energy-bin likelihoods of 14 years of Fermi-LAT data, their dark-matter TS
  profiles and their J-factor table), and the photon spectra per annihilation from PPPC 4 DM ID
  (Cirelli et al. 2011) and CosmiXs (Arina et al. 2024), via the gammapy-data mirror.
* **Direct detection**: LZ 2024 (HEPData), XENONnT 2025 (the collaboration's GitHub), earlier
  limits and the xenon neutrino fog compiled by C. O'Hare (GitHub).
* **Nuclear rates** for the neutron's fate: the PRIMAT compilation as shipped with PRyMordial.
* **Papers** (arXiv sources; only the named members are kept): TeX tables and vector figures.

Every file lands in data/ with its URL, size, SHA-256 and retrieval date in data/manifest.json.
GitHub files are pinned to a commit, so the bytes cannot change under us.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
import os
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import config as C  # noqa: E402

UA = "Mozilla/5.0 (X11; Linux x86_64) astrophysics-atlas-P4"
#: HEPData sits behind a bot filter that turns away browser-like strings it does not trust but
#: lets ordinary scripting clients through; we say what we are: a Python script.
UA_SCRIPT = "python-requests/2.32 astrophysics-atlas-P4"

GH_GAMMAPY = "https://raw.githubusercontent.com/gammapy/gammapy-data/43a8a09004975ce54342a440b077b7933485dabe"
GH_XENONNT = "https://raw.githubusercontent.com/XENONnT/wimp_data_release/0a5942061f0667c9b82d50aa6060b33f9c9b732e"
GH_OHARE = "https://raw.githubusercontent.com/cajohare/NeutrinoFog/0df3d0cd2322602dc147a66d94055f7f64879d80"
GH_PRYM = "https://raw.githubusercontent.com/vallima/PRyMordial/725d8a8db3ad5ea2630580d825c9d0d69ed76533"
SS = "https://member.ipmu.jp/satoshi.shirai"
FIGSHARE = "https://ndownloader.figshare.com/files"

#: (local path under data/raw, url, what it is, sibling cache to hard-link from, user agent)
FILES = [
    # --- equation of state ---------------------------------------------------------------
    ("eos/ss2018_standardmodel.dat", f"{SS}/standardmodel2018.dat",
     "Saikawa & Shirai 2018: T [GeV], g*rho, err, g*s, err", "68-thermal-history/data/raw/ss2018_standardmodel.dat", UA),
    ("eos/ss2020_eos2020.dat", f"{SS}/DM2020/eos2020.dat",
     "Saikawa & Shirai 2020: T [GeV], g_s, g_rho (the table behind their relic cross sections)", None, UA),
    ("eos/ss2020_eos2020_err1.dat", f"{SS}/DM2020/eos2020_err1.dat", "Saikawa & Shirai 2020: upper EoS curve", None, UA),
    ("eos/ss2020_eos2020_err2.dat", f"{SS}/DM2020/eos2020_err2.dat", "Saikawa & Shirai 2020: lower EoS curve", None, UA),
    ("eos/lms2015_standardmodel.dat", "http://www.laine.itp.unibe.ch/eos15/standardmodel.dat",
     "Laine & Meyer 2015 (Dec 2019 update): T/MeV, p/T^4, e/T^4, s/T^3, c/T^3, w, cs2, g_eff, h_eff, i_eff",
     "68-thermal-history/data/raw/lms2015_standardmodel.dat", UA),
    ("eos/escudero2020_sm_thermo.dat", "https://arxiv.org/src/2001.04466/anc/SM_Thermo.dat",
     "Escudero 2020 (JCAP 05, 048): t, T_gamma, T_gamma/T_nu, g*, g_s ... from 10 MeV to 1 keV",
     "68-thermal-history/data/raw/escudero2020_sm_thermo.dat", UA),
    # --- relic cross sections: answer keys ------------------------------------------------
    *[(f"relic/ss2020_{n}.dat", f"{SS}/DM2020/{n}.dat",
       f"Saikawa & Shirai 2020: m [GeV], sigma v [1e-26 cm^3/s], err (Omega), err (EoS); {n}", None, UA)
      for n in ("S_1", "S_2", "S_3", "S_4", "P_1", "P_2", "P_3", "P_4")],
    ("relic/sdb12_m-sv-num.txt", "https://arxiv.org/src/1204.3622v3/anc/m-sv-num",
     "Steigman, Dasgupta & Beacom 2012: m [GeV], <sigma v> [1e-26 cm^3/s] for Omega h^2 = 0.11", None, UA),
    # --- Fermi-LAT dwarfs ---------------------------------------------------------------
    ("fermi/mcdaniel2024_dSphs.tar.gz", f"{FIGSHARE}/42263853",
     "McDaniel et al. 2024 data release v2 (doi:10.6084/m9.figshare.24058650.v2): SEDs + TS profiles", None, UA),
    ("fermi/mcdaniel2024_dSphs.csv", f"{FIGSHARE}/44477168",
     "McDaniel et al. 2024: dwarf table (Table 1; J-factors corrected in v2)", None, UA),
    ("fermi/mcdaniel2024_basic_data_usage.ipynb", f"{FIGSHARE}/42768952",
     "McDaniel et al. 2024: the release's worked notebook (their flux and likelihood conventions)", None, UA),
    ("fermi/pppc4dmid_AtProduction_gammas.dat", f"{GH_GAMMAPY}/dark_matter_spectra/PPPC4DMID/AtProduction_gammas.dat",
     "PPPC 4 DM ID (Cirelli et al. 2011): dN/dlog10x of photons per annihilation, with EW corrections", None, UA),
    ("fermi/cosmixs_AtProduction-Gamma.dat", f"{GH_GAMMAPY}/dark_matter_spectra/cosmixs/AtProduction-Gamma.dat",
     "CosmiXs (Arina et al. 2024, JCAP 03, 035): photons per annihilation", None, UA),
    # --- direct detection -------------------------------------------------------------
    ("dd/lz2024_si.csv", "https://www.hepdata.net/download/table/ins2841863/SI%20cross%20section/1/csv",
     "LZ 2024 (PRL 135, 011802; HEPData ins2841863 v1, table 'SI cross section'): 90% CL limit + sensitivity", None, UA_SCRIPT),
    ("dd/xenonnt2025_si.csv", f"{GH_XENONNT}/SR0%2BSR1/xenonnt_2025_si_wimp.csv",
     "XENONnT SR0+SR1 2025 (arXiv:2502.18005): SI limit", None, UA),
    ("dd/ohare_Xe_SI_floor.txt", f"{GH_OHARE}/data/floors/Xe_SI.txt",
     "O'Hare 2021 (PRL 127, 251802): the xenon neutrino floor (edge of the fog)", None, UA),
    *[(f"dd/ohare_{n}", f"{GH_OHARE}/data/WIMPLimits/SI/{n}", f"O'Hare's compilation of SI limits: {n}", None, UA)
      for n in ("XENON1T-2018.txt", "LUX.txt", "PandaX.txt", "CDMSLite.txt", "CRESST.txt", "DarkSide.txt",
                "DAMA1.txt", "DAMA2.txt", "XENON1T.txt")],
    # --- the Higgs portal: width of a Standard-Model-like Higgs at any mass -----------------------
    ("higgs/yr4_Higgs_XSBR_update.xlsx",
     "https://twiki.cern.ch/twiki/pub/LHCPhysics/HiggsXSBR/Higgs_XSBR_YR4_update.xlsx",
     "LHC Higgs Working Group, CERN Yellow Report 4 (2016) workbook; sheet 'YR4 BSM Width': partial "
     "widths of a SM-like Higgs from 20 GeV to 1 TeV (without NLO EW corrections)", None, UA),
    # --- nuclear rates for the neutron's fate ---------------------------------------------
    *[(f"bbn/primat_{n}.txt", f"{GH_PRYM}/PRyMrates/nuclear/key_primat_rates/{n}.txt",
       f"PRIMAT rate (Pitrou et al. 2018) via PRyMordial: T9, N_A<sigma v>, uncertainty factor; {n}", None, UA)
      for n in ("npdg", "dpHe3g", "ddHe3n", "ddtp", "tpag", "tdan", "He3ntp", "He3dap")],
    ("bbn/prym_init.py", f"{GH_PRYM}/PRyM/PRyM_init.py",
     "PRyMordial's constants (reverse-rate coefficients we derive ourselves and compare with)", None, UA),
]

#: (arXiv id, members to keep, what it is). Only the named members are extracted.
ARXIV_SOURCES = [
    ("1204.3622", ["DMrelic-prd-rev1.tex"], "Steigman, Dasgupta & Beacom 2012 (PRD 86, 023506)"),
    ("2005.03544", ["main.tex", "anc/S_wave_Fig2.dat", "anc/P_wave_Fig3.dat"],
     "Saikawa & Shirai 2020 (JCAP 08, 011)"),
    ("2311.04982", ["main.tex", "dSph_table.tex", "result_shortcuts.tex"],
     "McDaniel et al. 2024 (PRD 109, 063024)"),
    ("1503.02641", ["main.tex", "p302_6yr_source_composite_limit_data_plot_bb.pdf",
                    "p302_6yr_source_composite_limit_data_plot_tautau.pdf"],
     "Ackermann et al. 2015 (PRL 115, 231301): the six-year limit figures (vector PDF)"),
    ("1306.4710", ["singlet-update12g.tex", "Fig1a.pdf", "Fig1b.pdf"],
     "Cline, Kainulainen, Scott & Weniger 2013 (PRD 88, 055025; v4 2015): scalar singlet formulas"),
    ("1411.4647", ["GCfit.tex"], "Calore, Cholis, McCabe & Weniger 2015 (PRD 91, 063003): GCE fits (Table)"),
]

#: Members of the McDaniel tarball we unpack (regular files only, names checked).
MCD_DIR = C.RAW_DIR / "fermi" / "mcdaniel2024"

#: PRyMordial (Burns, Tait & Valli 2024), the BBN code used to check the neutron calculation; the
#: same pinned commit as P5 (71-antimatter), whose zip is hard-linked when present.
PRYM_COMMIT = "725d8a8db3ad5ea2630580d825c9d0d69ed76533"
PRYM_URL = f"https://codeload.github.com/vallima/PRyMordial/zip/{PRYM_COMMIT}"
PRYM_ZIP = C.DATA_DIR / "external" / f"PRyMordial-{PRYM_COMMIT[:12]}.zip"
PRYM_DIR = C.DATA_DIR / "external" / f"PRyMordial-{PRYM_COMMIT}"
PRYM_SIBLING = f"71-antimatter/data/external/PRyMordial-{PRYM_COMMIT[:12]}.zip"


def fetch_prymordial(manifest: dict, today: str, force: bool = False, verbose: bool = True) -> None:
    import zipfile

    PRYM_ZIP.parent.mkdir(parents=True, exist_ok=True)
    key = f"external/{PRYM_ZIP.name}"
    if force or not PRYM_ZIP.exists():
        sib = C.ROOT / PRYM_SIBLING
        if sib.exists() and not force:
            if PRYM_ZIP.exists():
                PRYM_ZIP.unlink()
            os.link(sib, PRYM_ZIP)
            extra = {"hard_link_from": PRYM_SIBLING, "linked": today}
        else:
            if verbose:
                print(f"  fetching {PRYM_URL}")
            PRYM_ZIP.write_bytes(_get(PRYM_URL))
            extra = {"fetched": today}
        b = PRYM_ZIP.read_bytes()
        manifest[key] = {"url": PRYM_URL, "what": "PRyMordial source at a pinned commit "
                         "(github.com/vallima/PRyMordial)", "bytes": len(b), "sha256": _sha256(b), **extra}
    if force or not PRYM_DIR.exists():
        with zipfile.ZipFile(PRYM_ZIP) as z:
            for m in z.namelist():                          # refuse path traversal
                if m.startswith("/") or ".." in Path(m).parts:
                    raise RuntimeError(f"unsafe path in PRyMordial zip: {m}")
            z.extractall(PRYM_ZIP.parent)


def _get(url: str, ua: str = UA, *, tries: int = 4, timeout: float = 120.0) -> bytes:
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
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


def _looks_like_html(b: bytes) -> bool:
    head = b.lstrip()[:200].lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html")


def _unpack_mcdaniel(tar_path: Path, manifest: dict, today: str) -> None:
    """Extract SEDs/*.fits and TS_profiles/*.npy (and nothing else) from the release."""
    MCD_DIR.mkdir(parents=True, exist_ok=True)
    n = 0
    with tarfile.open(tar_path, mode="r:gz") as tf:
        for info in tf.getmembers():
            parts = Path(info.name).parts
            if not info.isfile() or len(parts) != 3 or parts[0] != "dSphs" or ".." in parts:
                continue
            sub, name = parts[1], parts[2]
            ok = (sub == "SEDs" and name.endswith("_sed.fits")) or (sub == "TS_profiles" and name.endswith(".npy"))
            if not ok or "/" in name:
                continue
            out = MCD_DIR / sub / name
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(tf.extractfile(info).read())
            n += 1
    manifest["raw/fermi/mcdaniel2024/"] = {"what": "unpacked from mcdaniel2024_dSphs.tar.gz",
                                          "files": n, "unpacked": today}


def fetch_all(force: bool = False, verbose: bool = True) -> dict:
    C.RAW_DIR.mkdir(parents=True, exist_ok=True)
    C.LIT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = _load_manifest()
    today = dt.date.today().isoformat()

    for rel, url, what, sibling, ua in FILES:
        path = C.RAW_DIR / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not force:
            continue
        src = C.ROOT / sibling if sibling else None
        if src is not None and src.exists() and not force:
            if path.exists():
                path.unlink()
            os.link(src, path)
            b = path.read_bytes()
            manifest[f"raw/{rel}"] = {"url": url, "what": what, "bytes": len(b), "sha256": _sha256(b),
                                      "hard_link_from": sibling, "linked": today}
            if verbose:
                print(f"  linked {sibling}")
            continue
        if verbose:
            print(f"  fetching {url}")
        b = _get(url, ua)
        if len(b) < 200 or _looks_like_html(b):
            raise RuntimeError(f"{url} did not return a data file ({len(b)} bytes)")
        path.write_bytes(b)
        manifest[f"raw/{rel}"] = {"url": url, "what": what, "bytes": len(b), "sha256": _sha256(b),
                                  "fetched": today}

    tar_path = C.RAW_DIR / "fermi/mcdaniel2024_dSphs.tar.gz"
    if force or not (MCD_DIR / "SEDs").exists():
        _unpack_mcdaniel(tar_path, manifest, today)

    fetch_prymordial(manifest, today, force=force, verbose=verbose)

    for arxiv_id, members, what in ARXIV_SOURCES:
        outdir = C.LIT_DIR / arxiv_id
        if all((outdir / m).exists() for m in members) and not force:
            continue
        url = f"https://arxiv.org/e-print/{arxiv_id}"
        if verbose:
            print(f"  fetching {url}")
        b = _get(url)
        with tarfile.open(fileobj=io.BytesIO(b), mode="r:gz") as tf:
            for m in members:
                info = tf.getmember(m)  # exact name only: nothing else is unpacked
                if not info.isfile():
                    raise RuntimeError(f"{arxiv_id}: {m} is not a regular file")
                data = tf.extractfile(info).read()
                out = outdir / m
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_bytes(data)
                manifest[f"literature/{arxiv_id}/{m}"] = {
                    "url": url, "what": what, "bytes": len(data), "sha256": _sha256(data),
                    "source_tarball_bytes": len(b), "fetched": today}

    C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    m = fetch_all(force=args.force)
    total = 0
    seen_tarballs = set()
    for k, v in m.items():
        if not isinstance(v, dict) or "hard_link_from" in v:
            continue
        if "source_tarball_bytes" in v:
            key = k.split("/")[1]
            if key not in seen_tarballs:
                seen_tarballs.add(key)
                total += v["source_tarball_bytes"]
        else:
            total += v.get("bytes", 0)
    print(f"  {len(m)} entries in the manifest, {total / 1e6:.1f} MB downloaded in all")
    return 0


if __name__ == "__main__":
    sys.exit(main())
