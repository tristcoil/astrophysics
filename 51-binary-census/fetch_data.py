"""Download and cache every input of the binary-star census.

    python fetch_data.py            # fetch whatever is missing (safe to re-run)

What it gets (all public, no accounts):

* Gaia DR3 (ESA archive): every source within 100 pc; all 443,205 two-body
  orbits, 338,215 acceleration solutions, 56,808 spectroscopic trend solutions
  and the binary_masses table, each joined to gaia_source; and the A/B stars
  within 500 pc, for wide companions of massive stars.
* VizieR: SB9 (spectroscopic binary orbits), Raghavan et al. 2010 (the census
  we check against), El-Badry et al. 2021's wide pairs, Hunt & Reffert 2023's
  open clusters (ages) and their members among Gaia's orbits (CDS xMatch).
* Mamajek's dwarf sequence and DEBCat (masses of eclipsing binaries).
* gaiamock (El-Badry et al. 2024) with its Gaia scanning-law files, plus a
  local build of GSL that its C kernel needs.

Every Gaia table is checked against a server-side COUNT(*), every download is
recorded in data/manifest.json with its query/URL, size and time.
"""

from __future__ import annotations

import fcntl
import hashlib
import io
import json
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import votable
from astropy.table import Table

import config as C
import gaia_tap as G

UA = {"User-Agent": "astrophysics-atlas-M2/1.0 (python urllib)"}

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _manifest() -> dict:
    if C.MANIFEST_PATH.exists():
        return json.loads(C.MANIFEST_PATH.read_text())
    return {}


def _record(key: str, **info) -> None:
    """Add one entry to data/manifest.json (locked: downloads may run in parallel)."""
    with open(C.MANIFEST_PATH.with_suffix(".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        m = _manifest()
        m[key] = {"retrieved_utc": G.now_utc(), **info}
        C.MANIFEST_PATH.write_text(json.dumps(m, indent=1, sort_keys=True, default=str) + "\n")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def http_get(url: str, *, timeout: float = 300.0, attempts: int = 4, data: bytes | None = None,
             headers: dict | None = None) -> bytes:
    last = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, data=data, headers={**UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # network trouble: retry with back-off
            last = exc
            time.sleep(5 * attempt)
    raise RuntimeError(f"GET {url[:120]} failed: {last}")


def download(url: str, path: Path, *, timeout: float = 600.0) -> Path:
    if path.exists():
        return path
    tmp = path.with_name(path.name + ".part")
    tmp.write_bytes(http_get(url, timeout=timeout))
    tmp.rename(path)
    return path


# ---------------------------------------------------------------------------
# Gaia DR3
# ---------------------------------------------------------------------------

#: gaia_source columns for the 100 pc sample. No quality cut is applied on the
#: archive side: binaries are exactly the stars a RUWE cut would throw away.
GAIA_100PC_COLUMNS = """
g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.parallax_over_error,
g.pmra, g.pmra_error, g.pmdec, g.pmdec_error, g.ruwe, g.astrometric_params_solved,
g.visibility_periods_used, g.astrometric_n_good_obs_al, g.astrometric_excess_noise,
g.ipd_gof_harmonic_amplitude, g.ipd_frac_multi_peak, g.phot_g_mean_mag, g.phot_bp_mean_mag,
g.phot_rp_mean_mag, g.bp_rp, g.phot_bp_rp_excess_factor, g.phot_g_mean_flux_over_error,
g.phot_bp_mean_flux_over_error, g.phot_rp_mean_flux_over_error, g.radial_velocity,
g.radial_velocity_error, g.rv_nb_transits, g.rv_method_used, g.rv_chisq_pvalue,
g.rv_renormalised_gof, g.rv_amplitude_robust, g.rv_template_teff, g.grvs_mag,
g.non_single_star, g.phot_variable_flag, g.teff_gspphot, s.fidelity_v2
""".strip()

NSS_ORBIT_COLUMNS = """
n.source_id, n.nss_solution_type, n.period, n.period_error, n.t_periastron, n.eccentricity,
n.eccentricity_error, n.arg_periastron, n.center_of_mass_velocity, n.center_of_mass_velocity_error,
n.semi_amplitude_primary, n.semi_amplitude_primary_error, n.semi_amplitude_secondary,
n.semi_amplitude_secondary_error, n.mass_ratio, n.mass_ratio_error,
n.a_thiele_innes, n.a_thiele_innes_error, n.b_thiele_innes, n.b_thiele_innes_error,
n.f_thiele_innes, n.f_thiele_innes_error, n.g_thiele_innes, n.g_thiele_innes_error,
n.inclination, n.inclination_error, n.parallax AS nss_parallax, n.parallax_error AS nss_parallax_error,
n.goodness_of_fit, n.efficiency, n.significance, n.flags, n.conf_spectro_period,
n.rv_n_good_obs_primary, n.rv_n_good_obs_secondary, n.astrometric_n_good_obs_al, n.g_luminosity_ratio,
g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.ruwe,
g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag, g.bp_rp, g.radial_velocity,
g.radial_velocity_error, g.rv_nb_transits, g.rv_template_teff, g.grvs_mag, g.teff_gspphot,
g.logg_gspphot, g.mh_gspphot, g.ebpminrp_gspphot
""".strip()

NSS_ACCEL_COLUMNS = """
n.source_id, n.nss_solution_type, n.significance, n.goodness_of_fit, n.accel_ra, n.accel_ra_error,
n.accel_dec, n.accel_dec_error, n.deriv_accel_ra, n.deriv_accel_ra_error, n.deriv_accel_dec,
n.deriv_accel_dec_error, n.parallax AS nss_parallax, n.parallax_error AS nss_parallax_error,
g.ra, g.dec, g.parallax, g.phot_g_mean_mag, g.bp_rp, g.ruwe
""".strip()

NSS_TREND_COLUMNS = """
n.source_id, n.nss_solution_type, n.mean_velocity, n.first_deriv_velocity,
n.first_deriv_velocity_error, n.second_deriv_velocity, n.second_deriv_velocity_error,
n.goodness_of_fit, g.ra, g.dec, g.parallax, g.phot_g_mean_mag, g.bp_rp, g.ruwe
""".strip()

GAIA_QUERIES = {
    # name: (SELECT ... body with {lo}/{hi}, COUNT(*) body with {lo}/{hi}, description)
    "gaia_100pc": (
        f"SELECT {GAIA_100PC_COLUMNS}\nFROM gaiadr3.gaia_source AS g\n"
        "LEFT OUTER JOIN external.gaiaedr3_spurious AS s ON s.source_id = g.source_id\n"
        f"WHERE g.parallax > {C.PARALLAX_MIN_100PC!r} AND g.source_id >= {{lo}} AND g.source_id < {{hi}}",
        f"SELECT COUNT(*) FROM gaiadr3.gaia_source AS g WHERE g.parallax > {C.PARALLAX_MIN_100PC!r} "
        "AND g.source_id >= {lo} AND g.source_id < {hi}",
        "Every Gaia DR3 source with parallax > 10 mas, with Rybizki et al. (2022) fidelity_v2",
    ),
    "nss_orbits": (
        f"SELECT {NSS_ORBIT_COLUMNS}\nFROM gaiadr3.nss_two_body_orbit AS n\n"
        "JOIN gaiadr3.gaia_source AS g ON g.source_id = n.source_id\n"
        "WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "SELECT COUNT(*) FROM gaiadr3.nss_two_body_orbit AS n WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "All Gaia DR3 two-body orbits (astrometric, spectroscopic, eclipsing, combined)",
    ),
    "nss_accel": (
        f"SELECT {NSS_ACCEL_COLUMNS}\nFROM gaiadr3.nss_acceleration_astro AS n\n"
        "JOIN gaiadr3.gaia_source AS g ON g.source_id = n.source_id\n"
        "WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "SELECT COUNT(*) FROM gaiadr3.nss_acceleration_astro AS n WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "All Gaia DR3 astrometric acceleration solutions (orbits longer than the 34-month baseline)",
    ),
    "nss_trend": (
        f"SELECT {NSS_TREND_COLUMNS}\nFROM gaiadr3.nss_non_linear_spectro AS n\n"
        "JOIN gaiadr3.gaia_source AS g ON g.source_id = n.source_id\n"
        "WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "SELECT COUNT(*) FROM gaiadr3.nss_non_linear_spectro AS n WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "All Gaia DR3 radial-velocity trend solutions",
    ),
    "nss_vim": (
        "SELECT n.source_id, n.nss_solution_type, n.significance FROM gaiadr3.nss_vim_fl AS n\n"
        "WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "SELECT COUNT(*) FROM gaiadr3.nss_vim_fl AS n WHERE n.source_id >= {lo} AND n.source_id < {hi}",
        "All Gaia DR3 variability-induced-mover solutions",
    ),
    "binary_masses": (
        "SELECT b.* FROM gaiadr3.binary_masses AS b WHERE b.source_id >= {lo} AND b.source_id < {hi}",
        "SELECT COUNT(*) FROM gaiadr3.binary_masses AS b WHERE b.source_id >= {lo} AND b.source_id < {hi}",
        "Gaia DR3 binary_masses: component masses where an orbit and an isochrone mass combine",
    ),
    "bright_500pc": (
        f"SELECT {GAIA_100PC_COLUMNS}\nFROM gaiadr3.gaia_source AS g\n"
        "LEFT OUTER JOIN external.gaiaedr3_spurious AS s ON s.source_id = g.source_id\n"
        f"WHERE g.parallax > {C.PARALLAX_MIN_BRIGHT!r} AND g.parallax_over_error > 10 AND g.bp_rp < 0.5 "
        f"AND g.phot_g_mean_mag + 5*LOG10(g.parallax) - 10 < {C.MG_MAX_BRIGHT!r} "
        "AND g.source_id >= {lo} AND g.source_id < {hi}",
        f"SELECT COUNT(*) FROM gaiadr3.gaia_source AS g WHERE g.parallax > {C.PARALLAX_MIN_BRIGHT!r} "
        f"AND g.parallax_over_error > 10 AND g.bp_rp < 0.5 "
        f"AND g.phot_g_mean_mag + 5*LOG10(g.parallax) - 10 < {C.MG_MAX_BRIGHT!r} "
        "AND g.source_id >= {lo} AND g.source_id < {hi}",
        "Blue (BP-RP < 0.5), luminous (M_G < 2) stars within 500 pc: B and A stars",
    ),
}


def gaia_path(name: str) -> Path:
    return C.GAIA_DIR / f"{name}.fits"


def fetch_gaia(names=None) -> None:
    for name, (template, count_template, description) in GAIA_QUERIES.items():
        if names and name not in names:
            continue
        path = gaia_path(name)
        if path.exists():
            continue
        print(f"Gaia DR3: {name} ({description})", flush=True)
        t0 = time.time()
        table, log = G.fetch_all_sky(name, template, quiet=True)
        expected = G.count_all_sky(count_template)
        if len(table) != expected:
            raise RuntimeError(f"{name}: got {len(table):,} rows, the archive counts {expected:,}")
        if len(np.unique(table["source_id"])) != len(table):
            print(f"  note: {name} has repeated source_ids", flush=True)
        G.write_with_provenance(table, path, {"query_template": template, "server_count": expected,
                                              "description": description})
        # the per-pixel chunks only exist to make an interrupted download resumable
        shutil.rmtree(C.CHUNK_DIR / name, ignore_errors=True)
        _record(f"gaia/{name}", rows=len(table), query=template, server_count=expected,
                bytes=path.stat().st_size, seconds=round(time.time() - t0, 1))
        print(f"  {len(table):,} rows (server count {expected:,}) in {time.time() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------------------
# VizieR
# ---------------------------------------------------------------------------

def vizier_table(source: str, *, columns: str = "**", constraints: dict | None = None) -> Table:
    """Fetch a VizieR table (or the rows matching ``constraints``) through ASU as a VOTable."""
    params = {"-source": source, "-out": columns, "-out.max": "unlimited", "-oc.form": "dec"}
    params.update(constraints or {})
    url = "https://vizier.cds.unistra.fr/viz-bin/votable?" + urllib.parse.urlencode(params)
    body = http_get(url, timeout=900)
    if b"OVERFLOW" in body[:20000]:
        raise RuntimeError(f"{source}: VizieR flagged the output as truncated")
    vot = votable.parse(io.BytesIO(body), verify="ignore")
    table = vot.get_first_table().to_table(use_names_over_ids=True)
    return table


VIZIER_JOBS = {
    # file stem: (VizieR table, columns, constraints, note)
    "sb9_main": (C.SB9_TABLES["main"], "**", None, "SB9 systems: positions, magnitudes, spectral types"),
    "sb9_orbits": (C.SB9_TABLES["orbits"], "**", None, "SB9 orbits: P, e, K1, K2, grade"),
    "sb9_alias": (C.SB9_TABLES["alias"], "**", None, "SB9 other designations"),
    "raghavan10_stars": (C.RAGHAVAN_TABLES["stars"], "**", None, "Raghavan+2010 sample (454 stars)"),
    "raghavan10_companions": (C.RAGHAVAN_TABLES["companions"], "**", None, "Raghavan+2010 Table 13"),
    "raghavan10_companion_masses": (C.RAGHAVAN_TABLES["companion_masses"], "**", None,
                                    "Raghavan+2010 Table 18"),
    "elbadry21_200pc": (
        C.ELBADRY21_TABLE,
        "Source1,Source2,RA_ICRS,DE_ICRS,RA2deg,DE2deg,Plx1,Plx2,e_Plx1,e_Plx2,pmRA1,pmRA2,pmDE1,pmDE2,"
        "e_pmRA1,e_pmRA2,e_pmDE1,e_pmDE2,RUWE1,RUWE2,Gmag1,Gmag2,BPmag1,BPmag2,RPmag1,RPmag2,"
        "IPDfmp1,IPDfmp2,IPDgofha1,IPDgofha2,theta,sepAU,R,BinType",
        {"Plx1": ">5"}, "El-Badry+2021 pairs whose star 1 has parallax > 5 mas (within 200 pc)"),
    "elbadry21_bright": (
        C.ELBADRY21_TABLE,
        "Source1,Source2,RA_ICRS,DE_ICRS,Plx1,Plx2,e_Plx1,e_Plx2,RUWE1,RUWE2,Gmag1,Gmag2,"
        "BPmag1,BPmag2,RPmag1,RPmag2,theta,sepAU,R,BinType",
        {"Plx1": ">2", "Gmag1": "<10.5"},
        "El-Badry+2021 pairs within 500 pc with a bright star 1 (G < 10.5)"),
    "hunt23_clusters": (C.HUNT23_CLUSTERS, "**", None, "Hunt & Reffert 2023 clusters"),
}


def fetch_vizier() -> None:
    for stem, (source, columns, constraints, note) in VIZIER_JOBS.items():
        path = C.VIZIER_DIR / f"{stem}.fits"
        if path.exists():
            continue
        print(f"VizieR: {stem} ({note})", flush=True)
        t0 = time.time()
        table = vizier_table(source, columns=columns, constraints=constraints)
        for name in table.colnames:  # FITS cannot store object columns
            if table[name].dtype == object:
                table[name] = np.array([str(v) for v in table[name]])
        table.meta.clear()
        table.write(path, format="fits", overwrite=True)
        _record(f"vizier/{stem}", source=source, constraints=constraints, rows=len(table),
                bytes=path.stat().st_size, note=note)
        print(f"  {len(table):,} rows in {time.time() - t0:.0f} s", flush=True)


def cds_xmatch(table: Table, cat2: str, radius_arcsec: float, ra: str = "ra", dec: str = "dec") -> Table:
    """Positional cross-match of a local table against a VizieR catalogue (CDS xMatch)."""
    buf = io.StringIO()
    table.write(buf, format="ascii.csv")
    boundary = "----m2binarycensus"
    fields = {"request": "xmatch", "distMaxArcsec": str(radius_arcsec), "RESPONSEFORMAT": "votable",
              "cat2": cat2, "colRA1": ra, "colDec1": dec}
    body = [f'--{boundary}\r\nContent-Disposition: form-data; name="{k}"\r\n\r\n{v}\r\n'
            for k, v in fields.items()]
    body.append(f'--{boundary}\r\nContent-Disposition: form-data; name="cat1"; filename="cat1.csv"\r\n'
                f"Content-Type: text/csv\r\n\r\n{buf.getvalue()}\r\n--{boundary}--\r\n")
    raw = http_get(C.CDS_XMATCH, data="".join(body).encode(), timeout=1800,
                   headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    out = votable.parse(io.BytesIO(raw), verify="ignore").get_first_table().to_table()
    for name in out.colnames:
        if out[name].dtype == object:
            out[name] = np.array([str(v) for v in out[name]])
    out.meta.clear()
    return out


def fetch_cluster_members() -> None:
    """Hunt & Reffert members among Gaia's orbits, matched with CDS xMatch.

    xMatch pairs positions; the analysis keeps only matches whose Gaia DR3
    source_id is identical, so the 1-arcsec radius just limits the candidates.
    """
    path = C.VIZIER_DIR / "hunt23_members_nss.fits"
    if path.exists():
        return
    orbits = Table.read(gaia_path("nss_orbits"))
    types = np.char.strip(np.asarray(orbits["nss_solution_type"]).astype(str))
    keep = np.isin(types, ["SB1", "SB1C", "SB2", "SB2C", "AstroSpectroSB1", "Orbital", "EclipsingSpectro"])
    if not keep.any():
        raise RuntimeError("no orbits selected for the cluster cross-match")
    sub = orbits[keep]["source_id", "ra", "dec"]
    print(f"CDS xMatch: {len(sub):,} Gaia orbits against Hunt & Reffert 2023 members", flush=True)
    table = cds_xmatch(sub, f"vizier:{C.HUNT23_MEMBERS}", 1.0)
    table.write(path, format="fits", overwrite=True)
    _record("vizier/hunt23_members_nss", service=C.CDS_XMATCH, cat2=C.HUNT23_MEMBERS, radius_arcsec=1,
            uploaded_rows=len(sub), rows=len(table))
    print(f"  {len(table):,} matches", flush=True)


def fetch_sb9_gaia() -> None:
    """Gaia DR3 counterparts of SB9 systems (CDS xMatch against VizieR I/355).

    SB9 positions are J2000 at various epochs; a 10-arcsec radius allows for
    proper motion. The analysis picks, among the candidates, the Gaia source
    whose G best matches the SB9 magnitude.
    """
    path = C.VIZIER_DIR / "sb9_gaia_xmatch.fits"
    if path.exists():
        return
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    main = Table.read(C.VIZIER_DIR / "sb9_main.fits")
    ok = np.array([len(str(r).strip()) > 0 for r in main["RAJ2000"]])
    main = main[ok]
    c = SkyCoord([f"{r.decode() if isinstance(r, bytes) else r} {d.decode() if isinstance(d, bytes) else d}"
                  for r, d in zip(main["RAJ2000"], main["DEJ2000"])], unit=(u.hourangle, u.deg))
    up = Table({"seq": np.asarray(main["Seq"]), "ra": c.ra.deg, "dec": c.dec.deg})
    print(f"CDS xMatch: {len(up):,} SB9 systems against Gaia DR3 (I/355/gaiadr3)", flush=True)
    table = cds_xmatch(up, "vizier:I/355/gaiadr3", 10.0)
    table.write(path, format="fits", overwrite=True)
    _record("vizier/sb9_gaia_xmatch", service=C.CDS_XMATCH, cat2="vizier:I/355/gaiadr3",
            radius_arcsec=10, uploaded_rows=len(up), rows=len(table))
    print(f"  {len(table):,} candidate matches", flush=True)


# ---------------------------------------------------------------------------
# Plain files
# ---------------------------------------------------------------------------

def fetch_files() -> None:
    for key, url, name in (("mamajek", C.MAMAJEK_URL, "mamajek_dwarf_sequence.txt"),
                           ("debcat", C.DEBCAT_URL, "debcat_debs.dat")):
        path = C.RAW_DIR / name
        if path.exists():
            continue
        print(f"Download: {name}", flush=True)
        download(url, path)
        _record(f"file/{key}", url=url, bytes=path.stat().st_size, sha256=_sha256(path))


# ---------------------------------------------------------------------------
# gaiamock + GSL
# ---------------------------------------------------------------------------

def fetch_gaiamock() -> None:
    d = C.GAIAMOCK_DIR
    for name in C.GAIAMOCK_FILES:
        path = d / name
        if not path.exists():
            download(f"{C.GAIAMOCK_REPO}/{C.GAIAMOCK_COMMIT}/{name}", path)
    (d / "COMMIT").write_text(C.GAIAMOCK_COMMIT + "\n")

    scans = d / "healpix_scans"
    if not (scans.exists() and len(list(scans.glob("healpix_64_*.fits"))) == 49152):
        zpath = d / "healpix_scans.zip"
        if not zpath.exists():
            print("gaiamock: downloading the scanning-law files (286 MB)", flush=True)
            download(C.GAIAMOCK_SCANS_URL, zpath, timeout=3600)
        if zpath.stat().st_size != C.GAIAMOCK_SCANS_BYTES or _sha256(zpath) != C.GAIAMOCK_SCANS_SHA256:
            raise RuntimeError("healpix_scans.zip does not match the expected size/SHA-256")
        print("gaiamock: unpacking 49,152 scan files", flush=True)
        with zipfile.ZipFile(zpath) as z:
            z.extractall(d)
        _record("gaiamock/healpix_scans", url=C.GAIAMOCK_SCANS_URL, bytes=C.GAIAMOCK_SCANS_BYTES,
                sha256=C.GAIAMOCK_SCANS_SHA256)

    so = d / "kepler_solve_astrometry.so"
    if so.exists():
        return
    gsl = C.BUILD_DIR / "gsl"
    if not (gsl / "lib" / "libgsl.a").exists():
        tarball = download(C.GSL_URL, C.BUILD_DIR / "gsl-2.8.tar.gz")
        if _sha256(tarball) != C.GSL_SHA256:
            raise RuntimeError("GSL tarball SHA-256 mismatch")
        print("Building GSL 2.8 (static, -fPIC) into data/build/gsl: a few minutes", flush=True)
        with tarfile.open(tarball) as tf:
            tf.extractall(C.BUILD_DIR, filter="data")
        src = C.BUILD_DIR / "gsl-2.8"
        subprocess.run(["./configure", f"--prefix={gsl}", "--disable-shared", "--enable-static",
                        "CFLAGS=-O2 -fPIC"], cwd=src, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["make", "-j6"], cwd=src, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["make", "install"], cwd=src, check=True, stdout=subprocess.DEVNULL)
    print("Compiling gaiamock's C kernel", flush=True)
    subprocess.run(["gcc", "-O2", "-shared", "-fPIC", "-o", str(so), str(d / "kepler_solve_astrometry.c"),
                    f"-I{gsl / 'include'}", str(gsl / "lib" / "libgsl.a"), str(gsl / "lib" / "libgslcblas.a"),
                    "-lm"], check=True)


def fetch_all() -> None:
    fetch_files()
    fetch_vizier()
    fetch_sb9_gaia()
    fetch_gaia()
    fetch_cluster_members()
    fetch_gaiamock()


if __name__ == "__main__":
    # e.g.  python fetch_data.py vizier sb9gaia   or   python fetch_data.py gaia:nss_orbits,nss_accel
    only = sys.argv[1:]
    if only:
        for step in only:
            if step.startswith("gaia:"):
                fetch_gaia(step[5:].split(","))
                continue
            {"files": fetch_files, "vizier": fetch_vizier, "gaia": fetch_gaia, "sb9gaia": fetch_sb9_gaia,
             "members": fetch_cluster_members, "gaiamock": fetch_gaiamock}[step]()
    else:
        fetch_all()
