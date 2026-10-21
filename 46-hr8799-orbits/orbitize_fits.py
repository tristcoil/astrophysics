"""Fits with the orbitize! package itself (Blunt et al. 2020).

    python orbitize_fits.py

* HR 8799 b, c, d, e one at a time with OFTI (Orbits For The Impatient,
  Blunt et al. 2017), the rejection sampler built for short orbital arcs, on
  Wang et al.'s data and priors. These single-planet posteriors are an
  independent check of our own joint "Unconstrained" fit (which differs only
  in sharing one parallax and star mass between the planets).
* beta Pic b with orbitize!'s parallel-tempered MCMC on the example data
  shipped with orbitize! (Nielsen et al. 2020 astrometry, 2003-2018). Its orbit
  is seen almost exactly edge-on, so the planet passes in front of and behind
  its star; the CRIRES radial velocity (Snellen et al. 2014: the planet was
  approaching us in December 2013) tells which passage is the one in front.
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd

import astrometry as A
import config as C

warnings.filterwarnings("ignore", category=RuntimeWarning)

OFTI_ORBITS = 3_000      # per planet: OFTI is slow on 1-mas GPI points (b: ~1.6 s per orbit)
#: beta Pic priors as in orbitize!'s beta Pic tutorials: Hipparcos parallax
#: (van Leeuwen 2007) and a 1.75 +- 0.05 Msun system mass.
BETAPIC_PLX = (51.44, 0.12)
BETAPIC_MASS = (1.75, 0.05)
BETAPIC_RV = (56643.0, -15.4, 1.7)          # MJD, planet RV relative to the star (km/s)
BETAPIC_MCMC = dict(num_temps=5, num_walkers=100, n_steps=300_000, burn_steps=2000, thin=10)


def _orbitize_table(df: pd.DataFrame, tag: str = "tmp"):
    """Our astrometry rows -> orbitize!'s input table (object 1)."""
    from astropy.table import Table
    from orbitize.read_input import read_file
    rows = []
    for _, r in df.iterrows():
        if r.kind == "radec":
            rows.append(dict(epoch=r.mjd, object=1, raoff=r.ra, raoff_err=r.ra_err, decoff=r.dec,
                             decoff_err=r.dec_err, radec_corr=r["corr"] if np.isfinite(r["corr"]) else 0.0))
        else:
            rows.append(dict(epoch=r.mjd, object=1, sep=r.sep, sep_err=r.sep_err, pa=r.pa, pa_err=r.pa_err))
    path = C.CACHE / f"orbitize_input_{tag}.csv"     # one file per fit: runs go side by side
    pd.DataFrame(rows).to_csv(path, index=False)
    return read_file(str(path))


def ofti_hr8799(planet: str, n_orbits: int = OFTI_ORBITS) -> np.ndarray:
    out = C.CACHE / f"ofti_W18_{planet}.npz"
    if out.exists():
        return np.load(out)["post"]
    import orbitize.driver
    df = A.dataset("W18")
    table = _orbitize_table(df[df.planet == planet], tag=f"W18_{planet}")
    drv = orbitize.driver.Driver(table, "OFTI", 1, C.W18_MSTAR[0], C.W18_PARALLAX[0],
                                 mass_err=C.W18_MSTAR[1], plx_err=C.W18_PARALLAX[1],
                                 system_kwargs={"tau_ref_epoch": C.TAU_REF_MJD})
    t0 = time.time()
    drv.sampler.run_sampler(n_orbits, num_samples=50_000, num_cores=1)
    post = drv.sampler.results.post
    labels = drv.system.labels
    np.savez_compressed(out, post=post, labels=np.array(labels))
    print(f"OFTI HR 8799 {planet}: {len(post)} orbits in {(time.time() - t0) / 60:.1f} min")
    return post


#: orbitize!'s MCMC for the planets whose arcs are too curved for OFTI (d, e):
#: OFTI accepts against the best likelihood seen so far, and with 1-mas GPI
#: points that running maximum stays far below the true peak, so it lets poor
#: orbits through (median ln L of its d samples -268, against -35 for b).
PLANET_MCMC = dict(num_temps=5, num_walkers=100, n_steps=300_000, burn_steps=2000, thin=10)


def mcmc_hr8799(planet: str, n_threads: int = 2) -> np.ndarray:
    out = C.CACHE / f"mcmc_W18_{planet}.npz"
    if out.exists():
        return np.load(out)["post"]
    import orbitize.driver
    df = A.dataset("W18")
    table = _orbitize_table(df[df.planet == planet], tag=f"W18mcmc_{planet}")
    kw = dict(PLANET_MCMC)
    n_steps, burn, thin = kw.pop("n_steps"), kw.pop("burn_steps"), kw.pop("thin")
    drv = orbitize.driver.Driver(table, "MCMC", 1, C.W18_MSTAR[0], C.W18_PARALLAX[0],
                                 mass_err=C.W18_MSTAR[1], plx_err=C.W18_PARALLAX[1],
                                 system_kwargs={"tau_ref_epoch": C.TAU_REF_MJD},
                                 mcmc_kwargs=dict(kw, num_threads=n_threads))
    t0 = time.time()
    _fork_start_method()
    drv.sampler.run_sampler(n_steps, burn_steps=burn, thin=thin)
    post = drv.sampler.results.post
    np.savez_compressed(out, post=post, labels=np.array(drv.system.labels))
    print(f"orbitize! MCMC HR 8799 {planet}: {len(post)} samples in {(time.time() - t0) / 60:.1f} min")
    return post


def _fork_start_method():
    """orbitize!'s MCMC always opens a multiprocessing pool; Python 3.14's default
    'forkserver' start method cannot re-import this module's callers, so use fork
    (as orbitize!'s own NestedSampler does)."""
    import multiprocessing as mp
    mp.set_start_method("fork", force=True)


def betapic(n_threads: int = 1):
    out = C.CACHE / "mcmc_betapic.npz"
    if out.exists():
        z = np.load(out)
        return z["post"], list(z["labels"])
    import orbitize.driver
    from orbitize.read_input import read_file
    table = read_file(str(C.RAW / C.ORBITIZE_BETAPIC))
    table = table[table["quant_type"] != "rv"]          # astrometry only; RV used for the branch below
    kw = dict(BETAPIC_MCMC)
    n_steps, burn, thin = kw.pop("n_steps"), kw.pop("burn_steps"), kw.pop("thin")
    drv = orbitize.driver.Driver(table, "MCMC", 1, BETAPIC_MASS[0], BETAPIC_PLX[0],
                                 mass_err=BETAPIC_MASS[1], plx_err=BETAPIC_PLX[1],
                                 mcmc_kwargs=dict(kw, num_threads=n_threads))
    t0 = time.time()
    _fork_start_method()
    drv.sampler.run_sampler(n_steps, burn_steps=burn, thin=thin)
    post = drv.sampler.results.post
    labels = list(drv.system.labels)
    np.savez_compressed(out, post=post, labels=np.array(labels),
                        tau_ref=drv.system.tau_ref_epoch)
    print(f"beta Pic b MCMC: {len(post)} samples in {(time.time() - t0) / 60:.1f} min")
    return post, labels


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2 and sys.argv[1] == "ofti":          # python orbitize_fits.py ofti b
        ofti_hr8799(sys.argv[2])
    elif len(sys.argv) > 2 and sys.argv[1] == "mcmc":          # python orbitize_fits.py mcmc e
        mcmc_hr8799(sys.argv[2], n_threads=int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    elif len(sys.argv) > 1 and sys.argv[1] == "betapic":
        betapic(n_threads=int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    else:
        for p in C.PLANETS:
            ofti_hr8799(p)
        betapic()
