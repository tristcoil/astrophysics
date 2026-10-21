"""Our own BAO measurement: the correlation function of DESI DR1 luminous red galaxies.

The baryon acoustic feature is an excess of galaxy pairs separated by the sound
horizon, ~100 h^-1 Mpc today. We measure it ourselves from DESI's first public
data release, in the same two redshift slices DESI used for its LRG results
(0.4 < z < 0.6 and 0.6 < z < 0.8):

1. **Positions.** Redshift -> comoving distance with a fiducial cosmology
   (Planck 2018, DESI's choice), then 3D Cartesian coordinates in h^-1 Mpc.
2. **Pair counts.** Galaxy-galaxy (DD), galaxy-random (DR) and random-random
   (RR) pairs in bins of separation s and of |mu|, the cosine of the angle to
   the line of sight, with each object's clustering weight. Our own C code
   (``paircount.c``, OpenMP, cell lists). RR uses the randoms split into
   chunks (Keihanen et al. 2019), which costs a fraction of the full count
   for a negligible loss of precision.
3. **The correlation function.** Landy & Szalay (1993):
   xi = (DD - 2 DR + RR) / RR with normalised counts, then the monopole
   xi_0(s) and quadrupole xi_2(s) by Legendre projection in mu.
4. **Error bars.** Delete-one jackknife over 100 sky regions; the Hartlap
   factor corrects the inverse covariance for the finite number of regions.
5. **The BAO fit.** xi_0(s) = B^2 xi_template(alpha s) + a0 + a1/s + a2/s^2
   over 50-150 h^-1 Mpc. The template is linear theory (CAMB) with the BAO
   wiggles damped by the galaxies' bulk motions (Sigma_nl = 8 h^-1 Mpc, the
   pre-reconstruction value). alpha is how much the ruler appears stretched
   relative to the fiducial cosmology, so D_V/r_d = alpha x (D_V/r_d)_fid.

Everything is cached under ``data/derived/bao/``.
"""

from __future__ import annotations

import ctypes
import json
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree

import background as bg
from config import DERIVED_DIR, PROJECT_DIR

OUT_DIR = DERIVED_DIR / "bao"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LIB_PATH = OUT_DIR / "libpaircount.so"

#: DESI DR1's fiducial cosmology (Planck 2018 base LCDM), as in DESI 2024 III.
FIDUCIAL = {"h": 0.6736, "ombh2": 0.02237, "omch2": 0.1200, "ns": 0.9649}
Z_BINS = [(0.4, 0.6), (0.6, 0.8)]
S_MAX, NS, NMU = 200.0, 100, 40     # 2 h^-1 Mpc bins, 40 bins in |mu|
N_REGIONS = 100                      # jackknife regions (split between the caps by area)
RR_SPLIT_SIZE = 2.0                  # each random chunk for RR ~ 2x the data
FIT_RANGE = (50.0, 150.0)
S_BIN = 4.0                          # final separation bins, h^-1 Mpc
SIGMA_NL = 8.0                       # pre-reconstruction BAO damping, h^-1 Mpc


# ---------------------------------------------------------------------------
# The pair counter
# ---------------------------------------------------------------------------


def load_library() -> ctypes.CDLL:
    src = PROJECT_DIR / "paircount.c"
    if not LIB_PATH.exists() or LIB_PATH.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC",
                        str(src), "-o", str(LIB_PATH), "-lm"], check=True)
    lib = ctypes.CDLL(str(LIB_PATH))
    d, i, L = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.c_long
    lib.count_pairs.argtypes = [d, d, i, d, L, d, d, i, d, L, ctypes.c_int, ctypes.c_double,
                                ctypes.c_int, ctypes.c_int, ctypes.c_int, d, d, d, d, ctypes.c_int]
    return lib


def _ptr(a, kind=ctypes.c_double):
    return a.ctypes.data_as(ctypes.POINTER(kind))


def pair_counts(cat1: dict, cat2: dict | None = None, *, nreg: int = N_REGIONS,
                s_max: float = S_MAX, ns: int = NS, nmu: int = NMU, threads: int = 0) -> dict:
    """Weighted pair counts; ``cat2=None`` means the auto-correlation of cat1."""
    lib = load_library()
    auto = cat2 is None
    c2 = cat1 if auto else cat2
    arrs = []
    for c in (cat1, c2):
        arrs.append((np.ascontiguousarray(c["pos"], float), np.ascontiguousarray(c["w"], float),
                     np.ascontiguousarray(c["reg"], np.int32), np.ascontiguousarray(c["z"], float)))
    total = np.zeros(ns * nmu)
    touch = np.zeros(nreg * ns * nmu)
    vsum, wsum = np.zeros(ns), np.zeros(ns)
    (p1, w1, r1, v1), (p2, w2, r2, v2) = arrs
    lib.count_pairs(_ptr(p1), _ptr(w1), _ptr(r1, ctypes.c_int), _ptr(v1), len(w1),
                    _ptr(p2), _ptr(w2), _ptr(r2, ctypes.c_int), _ptr(v2), len(w2),
                    int(auto), s_max, ns, nmu, nreg, _ptr(total), _ptr(touch), _ptr(vsum),
                    _ptr(wsum), threads)
    return {"total": total.reshape(ns, nmu), "touch": touch.reshape(nreg, ns, nmu),
            "vsum": vsum, "wsum": wsum}


def pair_norm(w: np.ndarray, reg: np.ndarray, auto: bool, w2=None, reg2=None, nreg=N_REGIONS):
    """Total pair weight (and its jackknife versions) that normalises a count."""
    W = np.bincount(reg, weights=w, minlength=nreg)
    if auto:
        W2 = np.bincount(reg, weights=w**2, minlength=nreg)
        tot = 0.5 * (W.sum() ** 2 - W2.sum())
        jk = 0.5 * ((W.sum() - W) ** 2 - (W2.sum() - W2))
        return tot, jk
    V = np.bincount(reg2, weights=w2, minlength=nreg)
    return W.sum() * V.sum(), (W.sum() - W) * (V.sum() - V)


# ---------------------------------------------------------------------------
# Catalogues
# ---------------------------------------------------------------------------


def fiducial_cosmo() -> bg.Cosmo:
    return bg.Cosmo(h=FIDUCIAL["h"], ombh2=FIDUCIAL["ombh2"], omch2=FIDUCIAL["omch2"],
                    w0=-1.0, wa=0.0)


def comoving_h(z: np.ndarray) -> np.ndarray:
    """Fiducial comoving distance in h^-1 Mpc, interpolated from a fine table."""
    zt = np.linspace(0.0, 1.5, 3001)
    grid = bg.DistanceGrid(zt[1:], z_max=1.5)
    dc = np.concatenate([[0.0], grid.comoving(fiducial_cosmo())[0]]) * FIDUCIAL["h"]
    return np.interp(z, zt, dc)


def cartesian(ra, dec, z) -> np.ndarray:
    r = comoving_h(z)
    ra, dec = np.radians(ra), np.radians(dec)
    return np.column_stack([r * np.cos(dec) * np.cos(ra), r * np.cos(dec) * np.sin(ra),
                            r * np.sin(dec)])


def _unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def raw_catalogues(cap: str) -> tuple[dict, dict]:
    """All LRGs and both random files of one cap, in a fixed row order."""
    import fetch_data as fd

    dat = fd.load_lss(cap, "data")
    rans = [fd.load_lss(cap, "random", i) for i in (0, 1)]
    return dat, {k: np.concatenate([r[k] for r in rans]) for k in rans[0]}


def load_catalogues(zmin: float, zmax: float, seed: int = 3, recon: bool = False) -> dict:
    """Data and randoms for both caps in one redshift slice, with jackknife regions.

    With ``recon`` each object also carries ``pos_rec``, its reconstructed
    position (``reconstruct_caps``); selection is always by observed redshift.
    """
    rng = np.random.default_rng(seed)
    out = {}
    caps = ("NGC", "SGC")
    raw = {}
    for cap in caps:
        dat, ran = raw_catalogues(cap)
        if recon:
            rec = reconstruct_caps()[cap]
            dat["pos_rec"], ran["pos_rec"] = rec["data"], rec["random"]
        raw[cap] = (dat, ran)
    n_ran = {cap: len(raw[cap][1]["ra"]) for cap in caps}
    n_reg = {"NGC": int(round(N_REGIONS * n_ran["NGC"] / sum(n_ran.values())))}
    n_reg["SGC"] = N_REGIONS - n_reg["NGC"]
    offset = {"NGC": 0, "SGC": n_reg["NGC"]}
    ratio = None
    for cap in caps:
        dat, ran = raw[cap]
        sel_d = (dat["z"] > zmin) & (dat["z"] < zmax)
        sel_r = (ran["z"] > zmin) & (ran["z"] < zmax)
        d = {k: v[sel_d] for k, v in dat.items()}
        r = {k: v[sel_r] for k, v in ran.items()}
        # jackknife regions: k-means on the sky positions of a random subsample
        sub = _unit(r["ra"], r["dec"])[rng.choice(len(r["ra"]), size=min(300000, len(r["ra"])),
                                                 replace=False)]
        cent, _ = kmeans2(sub, n_reg[cap], minit="++", seed=seed, iter=30)
        tree = cKDTree(cent)
        d["reg"] = tree.query(_unit(d["ra"], d["dec"]))[1].astype(np.int32) + offset[cap]
        r["reg"] = tree.query(_unit(r["ra"], r["dec"]))[1].astype(np.int32) + offset[cap]
        d["pos"] = cartesian(d["ra"], d["dec"], d["z"])
        r["pos"] = cartesian(r["ra"], r["dec"], r["z"])
        # same random-to-data weight ratio in both caps, so the caps combine cleanly
        if ratio is None:
            ratio = r["w"].sum() / d["w"].sum()
        r["w"] = r["w"] * ratio * d["w"].sum() / r["w"].sum()
        out[cap] = (d, r)
    return out


# ---------------------------------------------------------------------------
# Measuring xi(s, mu)
# ---------------------------------------------------------------------------


def measure(zmin: float, zmax: float, *, force: bool = False, threads: int = 0) -> dict:
    """DD, DR and RR in one redshift slice, kept per sky cap (cached).

    The two caps (north and south galactic) are thousands of Mpc apart, so no
    pair crosses between them: every count *and every normalisation* is within
    a cap, and the caps are combined afterwards (``combine``). RR is counted in
    chunks of the randoms; its normalisation is that of the chunks, and it is
    rescaled to the full random catalogue of its cap when combining.
    """
    path = OUT_DIR / f"paircounts_z{zmin:.1f}-{zmax:.1f}.npz"
    if path.exists() and not force:
        with np.load(path) as z:
            return {k: z[k] for k in z.files}
    t0 = time.monotonic()
    cats = load_catalogues(zmin, zmax)
    rng = np.random.default_rng(11)
    out = {"zmin": zmin, "zmax": zmax, "caps": np.array(list(cats))}
    for cap, (d, r) in cats.items():
        print(f"    {cap}: {len(d['w']):,} galaxies, {len(r['w']):,} randoms", flush=True)
        dd = pair_counts(d, threads=threads)
        dr = pair_counts(d, r, threads=threads)
        n_chunks = max(1, int(round(len(r["w"]) / (RR_SPLIT_SIZE * len(d["w"])))))
        lab = rng.integers(0, n_chunks, len(r["w"]))
        rr_tot, rr_touch, n_ch, n_ch_jk = 0.0, 0.0, 0.0, 0.0
        for k in range(n_chunks):
            sel = lab == k
            ch = {key: v[sel] for key, v in r.items()}
            rr = pair_counts(ch, threads=threads)
            rr_tot = rr_tot + rr["total"]
            rr_touch = rr_touch + rr["touch"]
            t, jk = pair_norm(ch["w"], ch["reg"], True)
            n_ch, n_ch_jk = n_ch + t, n_ch_jk + jk
        nDD = pair_norm(d["w"], d["reg"], True)
        nDR = pair_norm(d["w"], d["reg"], False, r["w"], r["reg"])
        nRR = pair_norm(r["w"], r["reg"], True)
        out.update({f"{cap}_DD": dd["total"], f"{cap}_DDt": dd["touch"],
                    f"{cap}_DR": dr["total"], f"{cap}_DRt": dr["touch"],
                    f"{cap}_RRc": rr_tot, f"{cap}_RRct": rr_touch,
                    f"{cap}_nDD": nDD[0], f"{cap}_nDDjk": nDD[1],
                    f"{cap}_nDR": nDR[0], f"{cap}_nDRjk": nDR[1],
                    f"{cap}_nRR": nRR[0], f"{cap}_nRRjk": nRR[1],
                    f"{cap}_nRRc": n_ch, f"{cap}_nRRcjk": n_ch_jk,
                    f"{cap}_vsum": dd["vsum"], f"{cap}_wsum": dd["wsum"],
                    f"{cap}_n_data": len(d["w"]), f"{cap}_n_random": len(r["w"]),
                    f"{cap}_n_chunks": n_chunks})
        print(f"      counted in {time.monotonic() - t0:.0f} s ({n_chunks} RR chunks)", flush=True)
    out["seconds"] = time.monotonic() - t0
    np.savez(path, **out)
    return out


RECON_PATH = OUT_DIR / "reconstructed_positions.npz"


def reconstruct_caps(force: bool = False, workers: int = 8) -> dict:
    """Reconstructed positions of every LRG and random point (all of 0.4 < z < 1.1),
    aligned with ``raw_catalogues`` row order; cached (float32)."""
    if RECON_PATH.exists() and not force:
        with np.load(RECON_PATH) as z:
            return {cap: {"data": z[f"{cap}_data"], "random": z[f"{cap}_random"]}
                    for cap in ("NGC", "SGC")}
    import reconstruction as R

    out, save = {}, {}
    for cap in ("NGC", "SGC"):
        dat, ran = raw_catalogues(cap)
        pd = cartesian(dat["ra"], dat["dec"], dat["z"])
        pr = cartesian(ran["ra"], ran["dec"], ran["z"])
        print(f"    reconstructing {cap}: {len(pd):,} galaxies, {len(pr):,} randoms", flush=True)
        sd, sr, info = R.reconstruct(pd, dat["w"], pr, ran["w"], workers=workers)
        out[cap] = {"data": (pd - sd).astype(np.float32), "random": (pr - sr).astype(np.float32)}
        save[f"{cap}_data"], save[f"{cap}_random"] = out[cap]["data"], out[cap]["random"]
        save[f"{cap}_rms_shift"] = info["rms_shift_data"]
    np.savez(RECON_PATH, **save)
    return out


def measure_recon(zmin: float, zmax: float, *, force: bool = False, threads: int = 0) -> dict:
    """Post-reconstruction counts: DD and DS (shifted data x shifted randoms) and SS,
    plus RR of the original randoms reused from ``measure`` (same chunks)."""
    path = OUT_DIR / f"paircounts_recon_z{zmin:.1f}-{zmax:.1f}.npz"
    if path.exists() and not force:
        with np.load(path) as z:
            return {k: z[k] for k in z.files}
    pre = measure(zmin, zmax)
    t0 = time.monotonic()
    cats = load_catalogues(zmin, zmax, recon=True)
    rng = np.random.default_rng(11)  # the same chunk labels as in ``measure``
    out = {"zmin": zmin, "zmax": zmax, "caps": np.array(list(cats))}
    for cap, (d, r) in cats.items():
        print(f"    {cap} (reconstructed): {len(d['w']):,} galaxies", flush=True)
        dS = dict(d, pos=d["pos_rec"].astype(float))
        rS = dict(r, pos=r["pos_rec"].astype(float))
        dd = pair_counts(dS, threads=threads)
        ds = pair_counts(dS, rS, threads=threads)
        n_chunks = int(pre[f"{cap}_n_chunks"])
        lab = rng.integers(0, n_chunks, len(r["w"]))
        ss_tot, ss_touch = 0.0, 0.0
        for k in range(n_chunks):
            sel = lab == k
            ch = {key: v[sel] for key, v in rS.items()}
            ss = pair_counts(ch, threads=threads)
            ss_tot, ss_touch = ss_tot + ss["total"], ss_touch + ss["touch"]
        keep = ("RRc", "RRct", "nDD", "nDDjk", "nDR", "nDRjk", "nRR", "nRRjk", "nRRc", "nRRcjk",
                "n_data", "n_random", "n_chunks")
        out.update({f"{cap}_{k}": pre[f"{cap}_{k}"] for k in keep})
        out.update({f"{cap}_DD": dd["total"], f"{cap}_DDt": dd["touch"],
                    f"{cap}_DR": ds["total"], f"{cap}_DRt": ds["touch"],
                    f"{cap}_SSc": ss_tot, f"{cap}_SSct": ss_touch,
                    f"{cap}_vsum": dd["vsum"], f"{cap}_wsum": dd["wsum"]})
        print(f"      counted in {time.monotonic() - t0:.0f} s", flush=True)
    out["seconds"] = time.monotonic() - t0
    np.savez(path, **out)
    return out


def combine(counts: dict, k: int | None = None) -> tuple:
    """Cap-combined DD, DR, RR (and SS) with normalisations; ``k`` drops jackknife region k.

    Returns the arguments of ``multipoles`` (SS and its normalisation are None
    before reconstruction).
    """
    DD = DR = RR = SS = 0.0
    nDD = nDR = nRR = nSS = 0.0
    recon = f"{counts['caps'][0]}_SSc" in counts
    for cap in counts["caps"]:
        cap = str(cap)
        g = lambda key: counts[f"{cap}_{key}"]
        if k is None:
            dd, dr, rrc = g("DD"), g("DR"), g("RRc")
            ndd, ndr, nrr, nrrc = g("nDD"), g("nDR"), g("nRR"), g("nRRc")
        else:
            dd, dr, rrc = g("DD") - g("DDt")[k], g("DR") - g("DRt")[k], g("RRc") - g("RRct")[k]
            ndd, ndr, nrr, nrrc = g("nDDjk")[k], g("nDRjk")[k], g("nRRjk")[k], g("nRRcjk")[k]
        DD, DR = DD + dd, DR + dr
        RR = RR + rrc * (nrr / nrrc)  # chunked RR scaled to the full randoms of this cap
        nDD, nDR, nRR = nDD + ndd, nDR + ndr, nRR + nrr
        if recon:  # shifted randoms carry the same weights, so the same full normalisation
            ssc = g("SSc") if k is None else g("SSc") - g("SSct")[k]
            nssc = g("nRRc") if k is None else g("nRRcjk")[k]
            SS = SS + ssc * (nrr / nssc)
            nSS = nSS + nrr
    if not recon:
        return DD, DR, RR, nDD, nDR, nRR, None, None
    return DD, DR, RR, nDD, nDR, nRR, SS, nSS


def multipoles(DD, DR, RR, nDD, nDR, nRR, SS=None, nSS=None,
               s_bin: float = S_BIN) -> tuple[np.ndarray, np.ndarray]:
    """Landy-Szalay xi(s, mu) -> monopole and quadrupole, in s bins of ``s_bin``.

    With SS (reconstruction) this is (DD - 2 DS + SS) / RR, where DR holds the
    data x shifted-random counts and RR the original randoms.
    """
    f = int(round(s_bin / (S_MAX / NS)))
    reb = lambda a: a.reshape(NS // f, f, NMU).sum(1)
    dd, dr, rr = reb(DD) / nDD, reb(DR) / nDR, reb(RR) / nRR
    ss = rr if SS is None else reb(SS) / nSS
    with np.errstate(invalid="ignore", divide="ignore"):
        xi = (dd - 2 * dr + ss) / rr
    mu = (np.arange(NMU) + 0.5) / NMU
    leg2 = 0.5 * (3 * mu**2 - 1)
    return xi.mean(axis=1), 5 * (xi * leg2).mean(axis=1)


def correlation_function(counts: dict, s_bin: float = S_BIN) -> dict:
    """xi_0, xi_2 and their delete-one jackknife covariances (separate and joint)."""
    xi0, xi2 = multipoles(*combine(counts), s_bin=s_bin)
    nreg = counts[f"{counts['caps'][0]}_DDt"].shape[0]
    jk = [multipoles(*combine(counts, k), s_bin=s_bin) for k in range(nreg)]
    jk0, jk2 = np.array([a for a, _ in jk]), np.array([b for _, b in jk])
    # delete-one jackknife: C = (N - 1)/N sum_k (x_k - mean)(x_k - mean)^T
    jcov = lambda a: (nreg - 1) / nreg * ((a - a.mean(0)).T @ (a - a.mean(0)))
    cov0, cov2 = jcov(jk0), jcov(jk2)
    cov02 = jcov(np.hstack([jk0, jk2]))
    # effective redshift: pair-weighted mean z of galaxy pairs at 80-120 h^-1 Mpc
    sl = slice(int(80 / (S_MAX / NS)), int(120 / (S_MAX / NS)))
    vs = sum(counts[f"{c}_vsum"][sl].sum() for c in counts["caps"])
    ws = sum(counts[f"{c}_wsum"][sl].sum() for c in counts["caps"])
    s = (np.arange(len(xi0)) + 0.5) * s_bin
    n_data = int(sum(counts[f"{c}_n_data"] for c in counts["caps"]))
    n_random = int(sum(counts[f"{c}_n_random"] for c in counts["caps"]))
    return {"s": s, "xi0": xi0, "xi2": xi2, "cov0": cov0, "cov2": cov2, "cov02": cov02,
            "n_regions": nreg, "z_eff": float(vs / ws), "n_data": n_data, "n_random": n_random}


# ---------------------------------------------------------------------------
# The BAO template and fit
# ---------------------------------------------------------------------------


def eh98_nowiggle(k, h, ombh2, omch2, ns):
    """Eisenstein & Hu (1998) no-wiggle transfer function shape, k in h/Mpc."""
    om, ob = ombh2 + omch2, ombh2
    fb = ob / om
    theta = 2.7255 / 2.7
    s = 44.5 * np.log(9.83 / om) / np.sqrt(1 + 10 * ob**0.75)
    a_g = 1 - 0.328 * np.log(431 * om) * fb + 0.38 * np.log(22.3 * om) * fb**2
    kmpc = k * h
    gamma = (om / h) * (a_g + (1 - a_g) / (1 + (0.43 * kmpc * s) ** 4))
    q = k * theta**2 / gamma
    L0 = np.log(2 * np.e + 1.8 * q)
    C0 = 14.2 + 731 / (1 + 62.5 * q)
    T = L0 / (L0 + C0 * q**2)
    return k**ns * T**2


def linear_power(k):
    """CAMB's linear matter power at z = 0 for the fiducial cosmology (h/Mpc units)."""
    path = OUT_DIR / "plin_fiducial.npz"
    if path.exists():
        with np.load(path) as z:
            return np.interp(np.log(k), np.log(z["k"]), z["p"])
    import camb

    pars = camb.set_params(H0=100 * FIDUCIAL["h"], ombh2=FIDUCIAL["ombh2"],
                           omch2=FIDUCIAL["omch2"], ns=FIDUCIAL["ns"], As=2.1e-9, mnu=0.06,
                           nnu=3.044, tau=0.0544, WantTransfer=True, kmax=50.0)
    pars.set_matter_power(redshifts=[0.0], kmax=50.0)
    res = camb.get_results(pars)
    kk, _, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=40.0, npoints=4000)
    np.savez(path, k=kk, p=pk[0])
    return np.interp(np.log(k), np.log(kk), pk[0])


def templates(s: np.ndarray, sigma_nl: float = SIGMA_NL):
    """xi_0 templates with and without the BAO feature, as functions of s."""
    k = np.logspace(-4, 1.3, 2**14)
    plin = linear_power(k)
    pnw = eh98_nowiggle(k, FIDUCIAL["h"], FIDUCIAL["ombh2"], FIDUCIAL["omch2"], FIDUCIAL["ns"])
    # rescale the EH shape to CAMB's broadband: smooth the ratio in ln k
    dlnk = np.log(k[1] / k[0])
    ratio = gaussian_filter1d(plin / pnw, sigma=0.25 / dlnk, mode="nearest")
    pnw = pnw * ratio
    pbao = (plin - pnw) * np.exp(-0.5 * k**2 * sigma_nl**2) + pnw
    damp = np.exp(-(k * 1.0) ** 2)  # 1 h^-1 Mpc smoothing so the transform converges

    def xi(p, ss):
        ks = np.outer(ss, k)
        j0 = np.sinc(ks / np.pi)
        return (j0 * (k**3 * p * damp)[None, :]).sum(1) * dlnk / (2 * np.pi**2)

    s_fine = np.linspace(20, 220, 801)
    return (interp1d(s_fine, xi(pbao, s_fine), kind="cubic"),
            interp1d(s_fine, xi(pnw, s_fine), kind="cubic"))


def fit_bao(cf: dict, *, sigma_nl: float = SIGMA_NL, fit_range=FIT_RANGE,
            alphas=np.linspace(0.80, 1.20, 801)) -> dict:
    """Profile chi^2(alpha) of xi_0 with B^2 and three broadband terms solved linearly."""
    s = cf["s"]
    m = (s >= fit_range[0]) & (s <= fit_range[1])
    ss, y = s[m], cf["xi0"][m]
    C = cf["cov0"][np.ix_(m, m)]
    nreg, nb = cf["n_regions"], m.sum()
    hartlap = (nreg - nb - 2) / (nreg - 1)
    icov = hartlap * np.linalg.inv(C)
    t_bao, t_nw = templates(s, sigma_nl)

    def scan(template):
        chi2 = np.empty(alphas.size)
        best = None
        for i, a in enumerate(alphas):
            X = np.column_stack([template(a * ss), np.ones_like(ss), 1 / ss, 1 / ss**2])
            F = X.T @ icov @ X
            beta = np.linalg.solve(F, X.T @ icov @ y)
            r = y - X @ beta
            chi2[i] = r @ icov @ r
            if best is None or chi2[i] < best[0]:
                best = (chi2[i], a, beta)
        return chi2, best

    chi2, best = scan(t_bao)
    chi2_nw, best_nw = scan(t_nw)
    i0 = int(np.argmin(chi2))
    inside = alphas[chi2 <= chi2[i0] + 1.0]
    s_model = np.linspace(fit_range[0], fit_range[1], 200)
    beta = best[2]
    model = (beta[0] * t_bao(best[1] * s_model) + beta[1] + beta[2] / s_model
             + beta[3] / s_model**2)
    beta_nw = best_nw[2]
    model_nw = (beta_nw[0] * t_nw(best_nw[1] * s_model) + beta_nw[1] + beta_nw[2] / s_model
                + beta_nw[3] / s_model**2)
    return {"alpha": float(best[1]), "alpha_lo": float(inside.min()),
            "alpha_hi": float(inside.max()), "chi2": float(chi2[i0]), "dof": int(nb - 5),
            "chi2_nowiggle": float(chi2_nw.min()),
            "detection_sigma": float(np.sqrt(max(chi2_nw.min() - chi2[i0], 0.0))),
            "alphas": alphas, "chi2_curve": chi2, "s_model": s_model, "model": model,
            "model_nowiggle": model_nw, "hartlap": hartlap, "sigma_nl": sigma_nl}


#: Anisotropic fit settings. Pre-reconstruction damping along / across the line of
#: sight and a small fingers-of-god scale (h^-1 Mpc), close to DESI DR1's LRG values.
SIGMA_PAR, SIGMA_PERP, SIGMA_FOG = 9.0, 4.5, 2.0
ANISO_S_BIN = 8.0                    # coarser bins keep the joint jackknife covariance usable
ANISO_RANGE = (48.0, 152.0)
BETAS = np.linspace(0.0, 1.0, 11)    # Kaiser beta = f/b, profiled


def aniso_templates(sigma_par=SIGMA_PAR, sigma_perp=SIGMA_PERP, sigma_fog=SIGMA_FOG):
    """Template multipoles xi_l(s) for l = 0, 2, 4, split by powers of beta.

    P(k, mu) = (1 + beta mu^2)^2 [P_nw + (P_lin - P_nw) exp(-k^2 (mu^2 S_par^2 +
    (1 - mu^2) S_perp^2) / 2)] / (1 + k^2 mu^2 S_fog^2 / 2)^2   (Kaiser; damped BAO; FoG)

    Returned as interpolators t[l][n] with xi_l = t[l][0] + beta t[l][1] + beta^2 t[l][2].
    """
    from scipy.special import eval_legendre, spherical_jn

    k = np.logspace(-4, 1.0, 6000)
    plin = linear_power(k)
    pnw = eh98_nowiggle(k, FIDUCIAL["h"], FIDUCIAL["ombh2"], FIDUCIAL["omch2"], FIDUCIAL["ns"])
    dlnk = np.log(k[1] / k[0])
    pnw = pnw * gaussian_filter1d(plin / pnw, sigma=0.25 / dlnk, mode="nearest")
    xg, wg = np.polynomial.legendre.leggauss(48)
    mu, wmu = 0.5 * (xg + 1), 0.5 * wg  # nodes on [0, 1]
    K, MU = np.meshgrid(k, mu, indexing="ij")
    damp = np.exp(-0.5 * K**2 * (MU**2 * sigma_par**2 + (1 - MU**2) * sigma_perp**2))
    base = (pnw[:, None] + (plin - pnw)[:, None] * damp) / (1 + 0.5 * (K * MU * sigma_fog) ** 2) ** 2
    s_fine = np.linspace(10, 260, 501)
    ks = np.outer(s_fine, k)
    soft = np.exp(-(k * 1.0) ** 2)
    out = {}
    for ell in (0, 2, 4):
        leg = eval_legendre(ell, mu)
        jl = spherical_jn(ell, ks)
        sign = {0: 1.0, 2: -1.0, 4: 1.0}[ell]
        out[ell] = []
        for power in (0, 1, 2):  # coefficient of beta^power in (1 + beta mu^2)^2
            coef = {0: 1.0, 1: 2.0, 2: 1.0}[power]
            pl = (2 * ell + 1) * ((base * coef * MU ** (2 * power) * leg[None, :]) @ wmu)
            xi = sign * (jl * (k**3 * pl * soft)[None, :]).sum(1) * dlnk / (2 * np.pi**2)
            out[ell].append(interp1d(s_fine, xi, kind="cubic", bounds_error=False,
                                     fill_value=0.0))
    return out


def ap_model(s, a_par, a_perp, beta, tmpl, nmu: int = 24):
    """Model xi_0, xi_2 in fiducial coordinates for dilations (alpha_par, alpha_perp).

    A separation (s, mu) measured with the fiducial cosmology is really
    s' = s sqrt(a_par^2 mu^2 + a_perp^2 (1 - mu^2)) at cosine mu' = mu a_par s / s'
    in the true one (Alcock & Paczynski 1979); project back onto Legendre multipoles.
    """
    xg, wg = np.polynomial.legendre.leggauss(nmu)
    mu, w = 0.5 * (xg + 1), 0.5 * wg
    fac = np.sqrt(a_par**2 * mu**2 + a_perp**2 * (1 - mu**2))
    sp = np.outer(s, fac)
    mup = (a_par * mu / fac)[None, :]
    xi = sum((t[0](sp) + beta * t[1](sp) + beta**2 * t[2](sp)) * lg
             for t, lg in ((tmpl[0], 1.0), (tmpl[2], 0.5 * (3 * mup**2 - 1)),
                           (tmpl[4], (35 * mup**4 - 30 * mup**2 + 3) / 8)))
    l2 = 0.5 * (3 * mu**2 - 1)
    return xi @ w, 5 * (xi * l2[None, :]) @ w


def fit_bao_aniso(cf: dict, *, sigmas=(SIGMA_PAR, SIGMA_PERP, SIGMA_FOG), fit_range=ANISO_RANGE,
                  grid=np.linspace(0.80, 1.20, 81)) -> dict:
    """Joint fit of xi_0 and xi_2 for alpha_par and alpha_perp.

    Linear parameters, solved exactly at every grid point: the amplitude B^2 and
    three broadband terms (a0 + a1/s + a2/s^2) for each multipole. beta is
    profiled over BETAS. alpha_par and alpha_perp get flat priors and are
    summarised from exp(-chi^2/2) on the grid.
    """
    s = cf["s"]
    m = (s >= fit_range[0]) & (s <= fit_range[1])
    ss = s[m]
    nb = m.sum()
    y = np.concatenate([cf["xi0"][m], cf["xi2"][m]])
    idx = np.concatenate([np.flatnonzero(m), len(s) + np.flatnonzero(m)])
    C = cf["cov02"][np.ix_(idx, idx)]
    nreg = cf["n_regions"]
    hartlap = (nreg - 2 * nb - 2) / (nreg - 1)
    icov = hartlap * np.linalg.inv(C)
    tmpl = aniso_templates(*sigmas)
    z = np.zeros_like(ss)
    poly0 = np.column_stack([np.ones_like(ss), 1 / ss, 1 / ss**2])
    P0 = np.vstack([poly0, np.zeros((nb, 3))])
    P2 = np.vstack([np.zeros((nb, 3)), poly0])
    chi2 = np.full((grid.size, grid.size), np.inf)
    best = (np.inf, None)
    for i, ap in enumerate(grid):
        for j, at in enumerate(grid):
            for beta in BETAS:
                x0, x2 = ap_model(ss, ap, at, beta, tmpl)
                X = np.column_stack([np.concatenate([x0, x2]), P0, P2])
                F = X.T @ icov @ X
                coef = np.linalg.solve(F, X.T @ icov @ y)
                r = y - X @ coef
                c = r @ icov @ r
                if c < chi2[i, j]:
                    chi2[i, j] = c
                    if c < best[0]:
                        best = (c, (ap, at, beta, coef))
    post = np.exp(-0.5 * (chi2 - chi2.min()))
    post /= post.sum()
    AP, AT = np.meshgrid(grid, grid, indexing="ij")
    mp, mt = (post * AP).sum(), (post * AT).sum()
    sp = np.sqrt((post * (AP - mp) ** 2).sum())
    st = np.sqrt((post * (AT - mt) ** 2).sum())
    rho = (post * (AP - mp) * (AT - mt)).sum() / (sp * st)
    # the isotropic dilation the two directions imply (DESI's D_V convention)
    aiso = AP ** (1 / 3) * AT ** (2 / 3)
    mi = (post * aiso).sum()
    si = np.sqrt((post * (aiso - mi) ** 2).sum())
    ap, at, beta, coef = best[1]
    s_model = np.linspace(fit_range[0], fit_range[1], 120)
    x0, x2 = ap_model(s_model, ap, at, beta, tmpl)
    pm = np.column_stack([np.ones_like(s_model), 1 / s_model, 1 / s_model**2])
    model0 = coef[0] * x0 + pm @ coef[1:4]
    model2 = coef[0] * x2 + pm @ coef[4:7]
    return {"alpha_par": float(mp), "alpha_par_err": float(sp), "alpha_perp": float(mt),
            "alpha_perp_err": float(st), "rho": float(rho), "alpha_iso": float(mi),
            "alpha_iso_err": float(si), "best_alpha_par": float(ap),
            "best_alpha_perp": float(at), "best_beta": float(beta), "chi2": float(best[0]),
            "dof": int(2 * nb - 2 - 1 - 7), "hartlap": hartlap, "chi2_grid": chi2,
            "grid": grid, "s_model": s_model, "model0": model0, "model2": model2,
            "edge": bool(post[[0, -1], :].sum() + post[:, [0, -1]].sum() > 0.01)}


def fiducial_distances(z: float) -> dict:
    """D_M/r_d, D_H/r_d and D_V/r_d of the fiducial cosmology at z."""
    from early_universe import EarlyUniverse

    c = fiducial_cosmo()
    grid = bg.DistanceGrid(np.array([z]), z_max=z)
    dm = grid.comoving(c)[0, 0]
    dh = bg.C_KM_S / bg.hubble(np.array([z]), c)[0, 0]
    rd = float(EarlyUniverse().rd(c)[0])
    return {"DM": dm / rd, "DH": dh / rd, "DV": float(np.cbrt(z * dm**2 * dh) / rd), "rd": rd}


def fiducial_dv_over_rd(z: float) -> tuple[float, float]:
    """(D_V / r_d) of the fiducial cosmology at z, and its r_d in Mpc."""
    from early_universe import EarlyUniverse

    c = fiducial_cosmo()
    grid = bg.DistanceGrid(np.array([z]), z_max=z)
    dm = grid.comoving(c)[0, 0]
    dh = bg.C_KM_S / bg.hubble(np.array([z]), c)[0, 0]
    rd = float(EarlyUniverse().rd(c)[0])
    return float(np.cbrt(z * dm**2 * dh) / rd), rd


#: After reconstruction the remaining smoothing is smaller (DESI DR1 LRG values, roughly).
POST_SIGMA_NL = 4.5
POST_SIGMAS = (6.0, 3.0, 2.0)


def analyse(counts: dict, sigma_iso: float, sigmas: tuple) -> tuple[dict, dict, dict]:
    """Isotropic and anisotropic fits of one set of pair counts."""
    cf = correlation_function(counts)
    cf8 = correlation_function(counts, s_bin=ANISO_S_BIN)
    iso = fit_bao(cf, sigma_nl=sigma_iso)
    an = fit_bao_aniso(cf8, sigmas=sigmas)
    fd = fiducial_distances(cf["z_eff"])
    sys_alpha = {f"sigma_nl={x:g}": fit_bao(cf, sigma_nl=x)["alpha"]
                 for x in (sigma_iso - 2.0, sigma_iso + 2.0)}
    sys_alpha["range 60-140"] = fit_bao(cf, sigma_nl=sigma_iso, fit_range=(60.0, 140.0))["alpha"]
    res = {"alpha": iso["alpha"], "alpha_lo": iso["alpha_lo"], "alpha_hi": iso["alpha_hi"],
           "dv_rd": iso["alpha"] * fd["DV"], "dv_rd_lo": iso["alpha_lo"] * fd["DV"],
           "dv_rd_hi": iso["alpha_hi"] * fd["DV"], "detection_sigma": iso["detection_sigma"],
           "chi2": iso["chi2"], "dof": iso["dof"], "sigma_nl": sigma_iso,
           "systematics_alpha": sys_alpha,
           "alpha_par": an["alpha_par"], "alpha_par_err": an["alpha_par_err"],
           "alpha_perp": an["alpha_perp"], "alpha_perp_err": an["alpha_perp_err"],
           "rho": an["rho"], "beta": an["best_beta"], "aniso_chi2": an["chi2"],
           "aniso_dof": an["dof"], "aniso_edge": an["edge"],
           "dm_rd": an["alpha_perp"] * fd["DM"], "dm_rd_err": an["alpha_perp_err"] * fd["DM"],
           "dh_rd": an["alpha_par"] * fd["DH"], "dh_rd_err": an["alpha_par_err"] * fd["DH"],
           "dv_rd_2d": an["alpha_iso"] * fd["DV"], "dv_rd_2d_err": an["alpha_iso_err"] * fd["DV"]}
    return res, {"cf": cf, "cf8": cf8, "iso": iso, "aniso": an}, fd


def run(*, force: bool = False, threads: int = 0, recon: bool = True) -> dict:
    """Measure and fit both redshift slices, before and (optionally) after reconstruction."""
    out = {}
    for zmin, zmax in Z_BINS:
        key = f"{zmin:.1f}-{zmax:.1f}"
        print(f"  DESI DR1 LRGs, {zmin} < z < {zmax}", flush=True)
        stages = [("pre", measure(zmin, zmax, force=force, threads=threads), SIGMA_NL,
                   (SIGMA_PAR, SIGMA_PERP, SIGMA_FOG))]
        if recon:
            stages.append(("post", measure_recon(zmin, zmax, force=force, threads=threads),
                           POST_SIGMA_NL, POST_SIGMAS))
        rec = {"zmin": zmin, "zmax": zmax}
        for stage, counts, siso, sig in stages:
            res, raw, fd = analyse(counts, siso, sig)
            rec.update({"z_eff": raw["cf"]["z_eff"], "n_data": raw["cf"]["n_data"],
                        "n_random": raw["cf"]["n_random"], "fiducial": fd})
            rec[stage] = res
            print(f"    {stage:4s}: alpha = {res['alpha']:.4f} (+{res['alpha_hi'] - res['alpha']:.4f}"
                  f"/-{res['alpha'] - res['alpha_lo']:.4f}), D_V/r_d = {res['dv_rd']:.3f}, "
                  f"BAO at {res['detection_sigma']:.1f} sigma; D_M/r_d = {res['dm_rd']:.2f} +- "
                  f"{res['dm_rd_err']:.2f}, D_H/r_d = {res['dh_rd']:.2f} +- {res['dh_rd_err']:.2f}, "
                  f"D_V/r_d (2D) = {res['dv_rd_2d']:.2f} +- {res['dv_rd_2d_err']:.2f}",
                  flush=True)
            cf, cf8, iso, an = raw["cf"], raw["cf8"], raw["iso"], raw["aniso"]
            np.savez(OUT_DIR / f"xi_{stage}_z{key}.npz", s=cf["s"], xi0=cf["xi0"], xi2=cf["xi2"],
                     cov0=cf["cov0"], cov2=cf["cov2"], s8=cf8["s"], xi0_8=cf8["xi0"],
                     xi2_8=cf8["xi2"], cov02_8=cf8["cov02"], s_model=iso["s_model"],
                     model=iso["model"], model_nowiggle=iso["model_nowiggle"],
                     alphas=iso["alphas"], chi2_curve=iso["chi2_curve"],
                     s_aniso=an["s_model"], model0_aniso=an["model0"], model2_aniso=an["model2"],
                     chi2_grid=an["chi2_grid"], grid=an["grid"])
        out[key] = rec
    (OUT_DIR / "bao_results.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Our own BAO measurement in DESI DR1 LRGs")
    ap.add_argument("--force", action="store_true", help="recount the pairs")
    ap.add_argument("--threads", type=int, default=0, help="OpenMP threads (0 = all cores)")
    ap.add_argument("--no-recon", action="store_true", help="skip reconstruction")
    args = ap.parse_args()
    run(force=args.force, threads=args.threads, recon=not args.no_recon)
