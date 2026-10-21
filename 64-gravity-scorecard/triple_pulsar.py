"""The triple pulsar PSR J0337+1715, simulated: how big is the signal of a falling-rate difference?

A millisecond pulsar and a white dwarf orbit each other every 1.63 days; that pair
orbits a second white dwarf every 327 days. Both members of the inner pair fall
toward the outer white dwarf with g ~ 1.7e-3 m/s^2. If the neutron star, whose
gravitational binding energy is ~10% of its mass, fell at a slightly different
rate (Delta = (m_G/m_I)_pulsar - (m_G/m_I)_WD != 0), the inner orbit would be
pushed toward (or away from) the outer star: a forced eccentricity that turns
with the outer orbit, seen in the pulse arrival times at frequency 2 f_I - f_O.

This is a labelled simulation with real initial conditions: the three bodies from
Ransom et al. (2014, Table 1: masses, osculating elements at MJD 55920; Newtonian
equations of motion, as they used), integrated with REBOUND/IAS15 with and without
the violation, in Voisin et al.'s (2020) form G_pi = G_po = G(1 + Delta), G_io = G.
The violation is switched on over three outer orbits, then the pulsar's line-of-sight
position (the Roemer delay) is followed for 2,988 days -- the span of the public
Nancay arrival times of Voisin et al. (2025). Everything an ordinary timing fit
would absorb (periods, phases, eccentricities, projected sizes, their slow drifts)
is fitted out; what is left at 2 f_I - f_O is the signature, which Archibald et
al. (2018) give as 30 ns per Delta = 1e-6.
"""

from __future__ import annotations

import json

import numpy as np
import rebound

import config as C

LIT = C.LITERATURE_DIR / "checks" / "pulsars.csv"
TIM = C.LIT_DATA_DIR / "papers" / "voisin-2025-zenodo" / "x" / "Data" / \
    "0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim"
CACHE = C.DERIVED_DIR / "triple_sim.npz"

K_GAUSS = 0.01720209895          # AU^(3/2) / (Msun^(1/2) day)
LT_S_AU = C.C_LIGHT / C.AU       # AU per light-second

#: Ransom et al. 2014, Table 1 (checked in literature/checks/pulsars.csv, rows tp_r14_*).
R14 = {
    "epoch_mjd": 55920.0,
    "m_p": 1.4378, "m_i": 0.19751, "m_o": 0.4101,
    "P_I": 1.629401788, "e_I": 6.9178e-4, "om_I": 97.6182, "tasc_I": 55920.407717436, "i_I": 39.254,
    "P_O": 327.257541, "e_O": 3.53561955e-2, "om_O": 95.619493, "tasc_O": 56233.935815, "i_O": 39.243,
    "dOmega": 2.7e-3, "x_I_lts": 1.21752844, "x_O_lts": 74.6727101,
}

RAMP_DAYS = 3 * 327.26
SPAN_DAYS = 2988.0
STEP_DAYS = 0.02
DELTA_SIM = 1e-5                 # exaggerated; the signal is linear in Delta (checked with +-)


def mean_anomaly_at(epoch, tasc, P, e, om_deg):
    """Mean anomaly (rad) at epoch for an orbit whose ascending node is passed at tasc."""
    nu = -np.radians(om_deg)                                  # true anomaly at the node (u = omega + nu = 0)
    E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
    M_asc = E - e * np.sin(E)
    return M_asc + 2 * np.pi / P * (epoch - tasc)


def build(delta: float) -> tuple[rebound.Simulation, callable]:
    p = R14
    sim = rebound.Simulation()
    sim.G = K_GAUSS**2
    sim.add(m=p["m_p"], name="pulsar")
    # inner white dwarf relative to the pulsar: the timing elements describe the pulsar's
    # orbit about the inner centre of mass, so the companion's periastron is 180 deg away
    sim.add(m=p["m_i"], P=p["P_I"], e=p["e_I"], omega=np.radians(p["om_I"] + 180), inc=np.radians(p["i_I"]),
            Omega=0.0, M=mean_anomaly_at(p["epoch_mjd"], p["tasc_I"], p["P_I"], p["e_I"], p["om_I"]), name="inner_wd")
    sim.add(m=p["m_o"], P=p["P_O"], e=p["e_O"], omega=np.radians(p["om_O"] + 180), inc=np.radians(p["i_O"]),
            Omega=np.radians(p["dOmega"]), M=mean_anomaly_at(p["epoch_mjd"], p["tasc_O"], p["P_O"], p["e_O"], p["om_O"]),
            name="outer_wd")
    sim.move_to_com()
    sim.integrator = "ias15"
    ps = sim.particles
    G = sim.G

    # Only the *difference* in how the pulsar and the inner white dwarf fall toward the outer
    # star is observable: a common change is the same as changing the masses, which a timing
    # fit absorbs (and which would otherwise drift the orbital phases by milliseconds). So the
    # pulsar's coupling to the outer star is raised by delta and the inner white dwarf's lowered
    # by delta*m_p/m_i: the inner pair's centre of mass then falls exactly as before, and the
    # differential acceleration is delta*(1 + m_p/m_i) times g_o. signature() rescales by that.
    eps = -delta * ps[0].m / ps[1].m

    def extra(simp):
        t = simp.contents.t
        s = 1.0 if t >= RAMP_DAYS else 0.5 * (1 - np.cos(np.pi * t / RAMP_DAYS))
        o = ps[2]
        for j, d in ((0, delta * s), (1, eps * s)):
            q = ps[j]
            dx, dy, dz = o.x - q.x, o.y - q.y, o.z - q.z
            r3 = (dx * dx + dy * dy + dz * dz) ** 1.5
            f = d * G / r3
            q.ax += f * o.m * dx
            q.ay += f * o.m * dy
            q.az += f * o.m * dz
            o.ax -= f * q.m * dx
            o.ay -= f * q.m * dy
            o.az -= f * q.m * dz

    if delta != 0.0:
        sim.additional_forces = extra
    return sim, extra


def run(delta: float):
    sim, _ = build(delta)
    times = np.arange(RAMP_DAYS, RAMP_DAYS + SPAN_DAYS, STEP_DAYS)
    z = np.empty(len(times))
    sim.integrate(RAMP_DAYS)
    for k, t in enumerate(times):
        sim.integrate(t)
        z[k] = sim.particles[0].z
    return times, z / LT_S_AU                    # line-of-sight position in light-seconds (Roemer delay)


def simulate(force: bool = False) -> dict:
    if CACHE.exists() and not force:
        z = np.load(CACHE)
        return {k: z[k] for k in z.files}
    out = {}
    for label, d in (("nominal", 0.0), ("plus", DELTA_SIM), ("minus", -DELTA_SIM)):
        t, z = run(d)
        out["t"], out[label] = t, z
        print(f"  triple run Delta={d:+.0e} done")
    np.savez(CACHE, **out)
    return out


def nuisance_basis(t_days, *, with_signature=True):
    """What a timing fit absorbs: Keplerian and three-body terms at the two orbital
    frequencies and their combinations, plus slow drifts; the signature last."""
    nI, nO = 2 * np.pi / R14["P_I"], 2 * np.pi / R14["P_O"]
    tt = (t_days - t_days.mean()) / SPAN_DAYS
    freqs = [nI, 2 * nI, 3 * nI, nO, 2 * nO, 3 * nO, nI + nO, nI - nO, 2 * nI + nO, nI + 2 * nO, nI - 2 * nO,
             2 * nI + 2 * nO, 2 * nI - 2 * nO]
    cols = [np.ones_like(tt), tt, tt**2]
    for w in freqs:
        cols += [np.cos(w * t_days), np.sin(w * t_days)]
    for w in (nI, 2 * nI, nO, 2 * nO):           # slow changes of period/phase/size
        cols += [tt * np.cos(w * t_days), tt * np.sin(w * t_days), tt**2 * np.cos(w * t_days), tt**2 * np.sin(w * t_days)]
    if with_signature:
        w = 2 * nI - nO
        cols += [np.cos(w * t_days), np.sin(w * t_days)]
    return np.column_stack(cols)


def signature(sim: dict | None = None) -> dict:
    sim = sim or simulate()
    t = sim["t"]
    d_plus = sim["plus"] - sim["nominal"]
    d_minus = sim["minus"] - sim["nominal"]
    odd = 0.5 * (d_plus - d_minus)                # linear part (seconds of Roemer delay)
    even = 0.5 * (d_plus + d_minus)
    A = nuisance_basis(t)
    x, *_ = np.linalg.lstsq(A, odd, rcond=None)
    resid = odd - A @ x
    amp = np.hypot(x[-2], x[-1])
    differential = DELTA_SIM * (1 + R14["m_p"] / R14["m_i"])     # see build(): pulsar minus inner WD
    per_1e6 = amp / differential * 1e-6
    scale = 1e-6 / differential
    return {"amplitude_ns_per_1e-6": per_1e6 * 1e9, "raw_peak_ns_per_1e-6": float(np.abs(odd).max() * scale * 1e9),
            "fit_residual_ns_per_1e-6": float(np.sqrt(np.mean(resid**2)) * scale * 1e9),
            "nonlinearity_ns_per_1e-6": float(np.abs(even).max() * scale * 1e9),
            "phase_rad": float(np.arctan2(x[-1], x[-2])),
            "published_A18_ns_per_1e-6": 30.0}


def load_toas(path=TIM):
    mjd, err = [], []
    for line in open(path):
        if line.startswith(("#", "C ", "FORMAT", "MODE")) or not line.strip():
            continue
        f = line.split()
        mjd.append(float(f[2]))
        err.append(float(f[3]))
    return np.array(mjd), np.array(err) * 1e-6


def white_noise_sensitivity(efac: float = 1.1177) -> dict:
    """Expected 1-sigma on Delta if the real Nancay TOAs (epochs, errors x EFAC) had white noise only."""
    mjd, err = load_toas()
    sig = signature()
    t = mjd - mjd.min() + RAMP_DAYS             # only phases matter; same basis
    A = nuisance_basis(t)
    w = 1 / (err * efac)
    cov = np.linalg.pinv((A * w[:, None]).T @ (A * w[:, None]))
    s_amp = np.sqrt(0.5 * (cov[-1, -1] + cov[-2, -2]))
    return {"n_toas": int(len(mjd)), "mjd_range": [float(mjd.min()), float(mjd.max())],
            "mean_err_us": float(err.mean() * 1e6), "median_err_us": float(np.median(err) * 1e6),
            "amp_sigma_ns": float(s_amp * 1e9),
            "delta_sigma_white": float(s_amp / (sig["amplitude_ns_per_1e-6"] * 1e-9) * 1e-6)}


if __name__ == "__main__":
    s = simulate()
    print(json.dumps(signature(s), indent=1))
    print(json.dumps(white_noise_sensitivity(), indent=1))
