"""Is our vacuum stable? The Higgs self-coupling up to the Planck scale.

With the couplings of ``rge.py`` matched to the measured masses and run up in
energy, this module finds

* where lambda turns negative (Lambda_lambda) and where the effective
  potential does (Lambda_I), its lowest value and the height of the barrier
  (Lambda_V = (max V_eff)^(1/4));
* the largest top mass for which the vacuum would be absolutely stable up to
  the Planck mass, and how many standard deviations the measured mass sits
  above it;
* the tunnelling rate of our vacuum at leading order -- a bubble of true
  vacuum of size R = 1/mu nucleates with probability ~ mu^4 exp(-S) per unit
  volume and time, S = 8 pi^2 / (3 |lambda(mu)|) (Coleman 1977; Isidori,
  Ridolfi & Strumia 2001) -- the probability that it has already happened in
  our past light cone, and lifetimes in the conventions of Buttazzo et al.
  (2013) and Andreassen, Frost & Schwartz (2018);
* the phase diagram in the (M_t, M_h) plane: stable, metastable, unstable;
* an independent check of matching and running with SMDR (Martin & Robertson
  2019), a C library that implements the state of the art.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import quad

import config as cfg
import rge

MU_MAX = 1e20                     # GeV, a little beyond the Planck mass
M_PL = cfg.M_PLANCK
SEC_PER_YEAR = 365.25 * 86400
SMDR_DIR = cfg.CACHE_DIR / "smdr" / "SMDR-1.3"

#: Planck 2018 (TT,TE,EE+lowE+lensing) for the light-cone volume and de Sitter future
H0_KMS_MPC = 67.36
OMEGA_L = 0.6847
OMEGA_R = 9.1e-5
MPC_KM = 3.0856775814913673e19


def h0_gev(h0=H0_KMS_MPC) -> float:
    return h0 / MPC_KM * cfg.HBAR_GEV_S


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def inputs(edition: str = "pdg2026", **override) -> rge.Inputs:
    """Measured inputs: PDG 2026 (default), PDG 2024 (the spec's values), or
    those of Buttazzo et al. (2013) and Andreassen et al. (2018), whose
    published numbers serve as benchmarks."""
    if edition == "pdg2026":
        d = dict(Mh=cfg.M_HIGGS[0], Mt=cfg.M_TOP[0], MW=cfg.M_W[0], MZ=cfg.M_Z[0],
                 alpha_s=cfg.ALPHA_S_MZ[0], GF=cfg.G_FERMI)
    elif edition == "pdg2024":
        d = dict(Mh=125.20, Mt=172.57, MW=80.3692, MZ=91.1880, alpha_s=0.1180, GF=1.1663788e-5)
    elif edition == "buttazzo2013":
        R = rge.REF
        d = dict(Mh=R["Mh"], Mt=R["Mt"], MW=R["MW"], MZ=R["MZ"], alpha_s=R["alpha_s"],
                 GF=1 / (np.sqrt(2) * R["V"] ** 2))
    elif edition == "afs2018":
        d = dict(Mh=125.09, Mt=173.1, MW=80.385, MZ=91.1876, alpha_s=0.1181, GF=1.1663787e-5)
    else:
        raise ValueError(edition)
    d.update(override)
    return rge.Inputs(**d)


# ---------------------------------------------------------------------------
# Trajectories and landmarks
# ---------------------------------------------------------------------------


COUPLINGS = ("lam", "yt", "g3", "g2", "gY", "yb", "ytau")


def smdr_offsets(edition: str = "pdg2026") -> dict | None:
    """SMDR's MS-bar couplings at M_t minus ours (Buttazzo et al.'s formulas),
    at the measured point. Added to our matching, they carry SMDR's state of
    the art (four-loop QCD in the top pole-mass relation, full electroweak
    two-loop matching, g3 from alpha_s at higher order) across a grid of
    masses, where calling SMDR (~30 s a point) would be too slow."""
    s = smdr_fit(inputs(edition), edition)
    if s is None:
        return None
    c = rge.match(inputs(edition))
    return {k: float(s[k] - c[k]) for k in COUPLINGS}


def match(inp: rge.Inputs, *, order: str = "NNLO", offsets: dict | None = None) -> dict:
    c = rge.match(inp, order)
    if offsets:
        c = dict(c)
        for k, v in offsets.items():
            c[k] = c[k] + v
    return c


def trajectory(inp: rge.Inputs, *, loops: int = 3, order: str = "NNLO",
               mu_max: float = MU_MAX, couplings: dict | None = None,
               offsets: dict | None = None) -> dict:
    c = couplings or match(inp, order=order, offsets=offsets)
    mu, X = rge.run(rge.state(c), c["mu"], mu_max, loops=loops)
    return {"mu": mu, "X": X, "start": c}


def _crossing(mu, y):
    """First scale where y changes sign from + to - (log-interpolated)."""
    idx = np.nonzero((y[:-1] > 0) & (y[1:] <= 0))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    t = y[i] / (y[i] - y[i + 1])
    return float(np.exp(np.log(mu[i]) + t * (np.log(mu[i + 1]) - np.log(mu[i]))))


def landmarks(tr: dict, *, eff_loops: int = 2) -> dict:
    mu, X = tr["mu"], tr["X"]
    lam = X[:, 6]
    big = mu > 1e4                              # h >> v: where lambda_eff applies
    leff = np.full(len(mu), np.nan)
    leff[big] = rge.lambda_eff(X[big].T, loops=eff_loops)
    below = mu <= M_PL
    i_min = int(np.argmin(np.where(below, lam, np.inf)))
    j = big & below
    V = np.where(j, leff * mu**4 / 4, -np.inf)
    i_v = int(np.argmax(np.where(leff > 0, V, -np.inf)))
    def at(m):
        return {n: float(np.interp(np.log(m), np.log(mu), X[:, k]))
                for k, n in enumerate(rge.NAMES)}
    pl = at(M_PL)
    return {"Lambda_lambda": _crossing(mu, lam),
            "Lambda_I": _crossing(mu[big], leff[big]),
            "Lambda_V": float(V[i_v] ** 0.25) if np.isfinite(V[i_v]) and V[i_v] > 0 else np.nan,
            "lambda_min": float(lam[i_min]), "mu_lambda_min": float(mu[i_min]),
            "lambda_eff_min": float(np.nanmin(np.where(j, leff, np.nan))),
            "lambda_MPl": pl["lam"], "yt_MPl": float(np.sqrt(pl["ytsq"])),
            "g1_MPl": float(np.sqrt(pl["g1sq"])), "g2_MPl": float(np.sqrt(pl["g2sq"])),
            "g3_MPl": float(np.sqrt(pl["g3sq"]))}


# ---------------------------------------------------------------------------
# Absolute stability bound on the top mass
# ---------------------------------------------------------------------------


def min_lambda_eff(inp: rge.Inputs, *, eff_loops: int = 2) -> float:
    tr = trajectory(inp, mu_max=M_PL)
    mu, X = tr["mu"], tr["X"]
    big = mu > 1e4
    return float(np.min(rge.lambda_eff(X[big].T, loops=eff_loops)))


def min_lambda(inp: rge.Inputs) -> float:
    tr = trajectory(inp, mu_max=M_PL)
    return float(np.min(tr["X"][:, 6]))


def run_many(inps: list[rge.Inputs], *, mu_max: float = M_PL, step: float = 0.05,
             offsets: dict | None = None):
    """Match each input set at its own M_t, carry all to a common starting
    scale, and run them together (vectorised). Returns mu, X (n, 8, N)."""
    starts = [match(i, offsets=offsets) for i in inps]
    mu0 = min(c["mu"] for c in starts)
    cols = []
    for c in starts:
        X = rge.state(c)
        if c["mu"] != mu0:
            X = rge.run(X, c["mu"], mu0, step=0.02)[1][-1]
        cols.append(X)
    return rge.run(np.stack(cols, axis=1), mu0, mu_max, step=step)


def criticality(inps: list[rge.Inputs], *, offsets: dict | None = None) -> dict:
    """min over 10^4 GeV < h < M_Pl of lambda_eff (two-loop) and of lambda."""
    mu, X = run_many(inps, offsets=offsets)
    big = mu > 1e4
    leff = np.stack([rge.lambda_eff(X[i], loops=2) for i in np.nonzero(big)[0]])
    return {"lambda_eff_min": leff.min(axis=0), "lambda_min": X[:, 6, :].min(axis=0)}


def stability_bound_mt(inp: rge.Inputs, *, criterion: str = "lambda_eff",
                       grid=np.arange(168.0, 176.01, 0.25), offsets: dict | None = None) -> float:
    """Top mass at which the vacuum becomes exactly critical (min over h < M_Pl
    of lambda_eff, or of lambda, equals zero), other inputs held fixed."""
    crit = criticality([_with(inp, Mt=float(t)) for t in grid], offsets=offsets)
    y = crit["lambda_eff_min" if criterion == "lambda_eff" else "lambda_min"]
    idx = np.nonzero((y[:-1] > 0) & (y[1:] <= 0))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    # cubic interpolation through four neighbours for the root
    j = slice(max(i - 1, 0), min(i + 3, len(grid)))
    p = np.polyfit(grid[j] - grid[i], y[j], min(3, len(grid[j]) - 1))
    roots = [r.real for r in np.roots(p) if abs(r.imag) < 1e-9 and -0.5 <= r.real <= 0.75]
    return float(grid[i] + roots[0]) if roots else float(np.interp(0, -y[i:i + 2], grid[i:i + 2]))


def _with(inp: rge.Inputs, **kw) -> rge.Inputs:
    d = dict(Mh=inp.Mh, Mt=inp.Mt, MW=inp.MW, MZ=inp.MZ, alpha_s=inp.alpha_s, GF=inp.GF,
             mb_mb=inp.mb_mb, Mtau=inp.Mtau)
    d.update(kw)
    return rge.Inputs(**d)


# ---------------------------------------------------------------------------
# Vacuum decay at leading order
# ---------------------------------------------------------------------------


def light_cone_volume(h0=None, om_l=OMEGA_L, om_r=OMEGA_R) -> float:
    """H0^4 x the space-time volume of our past light cone,
    (4 pi / 3) int a^4 (eta0 - eta)^3 d eta, for flat LambdaCDM
    (Buttazzo et al. eq. 6.3 quote ~0.15)."""
    om_m = 1 - om_l - om_r

    def E(a):
        return np.sqrt(om_r / a**4 + om_m / a**3 + om_l)

    def eta(a):                             # conformal time in units of 1/H0
        return quad(lambda x: 1 / (x * x * E(x)), 0, a, limit=200)[0]

    eta0 = eta(1.0)
    # d(eta) = da / (a^2 E): integrand a^4 (eta0 - eta)^3 / (a^2 E)
    val = quad(lambda a: a * a * (eta0 - eta(a)) ** 3 / E(a), 0, 1, limit=200)[0]
    return 4 * np.pi / 3 * val


def decay(tr: dict, *, light_cone: float | None = None) -> dict:
    """Leading-order nucleation rate per volume, Gamma/V = max_mu mu^4 e^{-S(mu)}
    over mu <= M_Pl, and what it implies."""
    mu, lam = tr["mu"], tr["X"][:, 6]
    ok = (mu <= M_PL) & (lam < 0)
    if not np.any(ok):
        return {"stable": True}
    log10_rate = np.where(ok, 4 * np.log10(mu) - rge.bounce_action(lam) / np.log(10), -np.inf)
    i = int(np.argmax(log10_rate))
    lr = float(log10_rate[i])                       # log10(Gamma/V / GeV^4)
    S = float(rge.bounce_action(lam[i]))
    H0 = h0_gev()
    lc = light_cone if light_cone is not None else light_cone_volume()
    log10_p0 = lr + np.log10(lc) - 4 * np.log10(H0)
    # Andreassen et al.: tau = (Gamma/V)^(-1/4)
    log10_tau_afs = -lr / 4 + np.log10(cfg.HBAR_GEV_S / SEC_PER_YEAR)
    # Buttazzo et al. eq. (6.6): future dominated by vacuum energy, H = H0 sqrt(Omega_L)
    Hds = H0 * np.sqrt(OMEGA_L)
    log10_tau_ds = (np.log10(3 / (4 * np.pi)) + 3 * np.log10(Hds) - lr
                    + np.log10(cfg.HBAR_GEV_S / SEC_PER_YEAR))
    return {"stable": False, "log10_rate_GeV4": lr, "bounce_scale": float(mu[i]),
            "lambda_at_bounce": float(lam[i]), "action": S,
            "light_cone_H0^4_volume": lc, "log10_p_decayed_in_past_light_cone": float(log10_p0),
            "log10_tau_years_afs": float(log10_tau_afs),
            "log10_tau_years_de_sitter": float(log10_tau_ds)}


# ---------------------------------------------------------------------------
# Phase diagram (vectorised over a grid of top and Higgs masses)
# ---------------------------------------------------------------------------


def phase_grid(mt: np.ndarray, mh: np.ndarray, base: rge.Inputs, *, step: float = 0.1,
               offsets: dict | None = None) -> dict:
    MT, MH = np.meshgrid(mt, mh)
    states = []
    for t, h in zip(MT.ravel(), MH.ravel()):
        states.append(rge.state(match(_with(base, Mt=float(t), Mh=float(h)), offsets=offsets)))
    X0 = np.stack(states, axis=1)                         # (8, N)
    # carry every column from its own M_t to a common starting scale, then
    # run all of them together
    mu0 = float(mt.min())
    Xc = np.empty_like(X0)
    for k, t in enumerate(MT.ravel()):
        _, Xk = rge.run(X0[:, k], float(t), mu0, step=0.02)
        Xc[:, k] = Xk[-1]
    mu, X = rge.run(Xc, mu0, M_PL, step=step)              # (n, 8, N)
    lam = X[:, 6, :]
    big = mu > 1e4
    leff = np.stack([rge.lambda_eff(X[i], loops=2) for i in np.nonzero(big)[0]])
    lmin = lam.min(axis=0)
    leffmin = leff.min(axis=0)
    # leading-order decay probability in the past light cone
    lc = light_cone_volume()
    H0 = h0_gev()
    with np.errstate(divide="ignore", over="ignore"):
        lr = np.where(lam < 0, 4 * np.log10(mu)[:, None] - rge.bounce_action(lam) / np.log(10),
                      -np.inf).max(axis=0)
    log10_p0 = lr + np.log10(lc) - 4 * np.log10(H0)
    lam_I = np.full(lam.shape[1], np.nan)
    for k in range(lam.shape[1]):
        lam_I[k] = _crossing(mu[big], leff[:, k])
    # non-perturbative: lambda above 1 anywhere below M_Pl, or a run that blew up
    # (finite steps can jump across a Landau pole and land on finite garbage)
    finite = np.all(np.isfinite(X[:, :7, :]), axis=(0, 1))
    lmax = np.where(finite, np.nanmax(np.where(np.isfinite(lam), lam, np.inf), axis=0), np.inf)
    return {"mt": mt, "mh": mh, "lambda_min": lmin.reshape(MT.shape),
            "lambda_max": lmax.reshape(MT.shape),
            "lambda_eff_min": leffmin.reshape(MT.shape),
            "log10_p0": log10_p0.reshape(MT.shape),
            "log10_Lambda_I": np.log10(lam_I).reshape(MT.shape),
            "lambda_MPl": lam[-1].reshape(MT.shape)}


# ---------------------------------------------------------------------------
# SMDR (Martin & Robertson 2019): independent matching and running
# ---------------------------------------------------------------------------

SMDR_INPUT_TEMPLATE = """SMDR_Mt_pole = {Mt};
SMDR_Mh_pole = {Mh};
SMDR_MZ_PDG = {MZ};
SMDR_alphaS_5_MZ = {alpha_s};
SMDR_alpha = 1/137.035999178;
SMDR_Delta_alpha_had_5_MZ_in = 0.02783;
SMDR_GFermi = {GF};
SMDR_mbmb = 4.183;
SMDR_mcmc = 1.273;
SMDR_ms_2GeV = 0.0935;
SMDR_md_2GeV = 0.00470;
SMDR_mu_2GeV = 0.00216;
SMDR_Mtau_pole = 1.77693;
SMDR_Mmuon_pole = 0.1056583755;
SMDR_Melectron_pole = 0.00051099895;
"""


def smdr_available() -> bool:
    return (SMDR_DIR / "calc_fit").exists() and (SMDR_DIR / "fig_RGrun_vs_Q").exists()


def smdr_fit(inp: rge.Inputs, tag: str) -> dict | None:
    """MS-bar parameters at Q = M_t from SMDR's calc_fit (its own G_F
    convention: the PDG value without the (3/10) M_mu^2/M_W^2 term)."""
    if not smdr_available():
        return None
    work = cfg.CACHE_DIR / "smdr_runs"
    work.mkdir(exist_ok=True)
    out_model = work / f"{tag}_model.dat"
    cache = work / f"{tag}_fit.json"
    if cache.exists():
        return json.loads(cache.read_text())
    gf = inp.GF * (1 - 0.3 * (0.1056583755 / inp.MW) ** 2)
    infile = work / f"{tag}_inputs.dat"
    infile.write_text(SMDR_INPUT_TEMPLATE.format(Mt=inp.Mt, Mh=inp.Mh, MZ=inp.MZ,
                                                 alpha_s=inp.alpha_s, GF=f"{gf:.10e}"))
    # SMDR keeps file names in 50-character buffers: pass short relative paths
    rel_in = Path("../../smdr_runs") / infile.name
    rel_out = Path("../../smdr_runs") / out_model.name
    res = subprocess.run(["./calc_fit", "-i", str(rel_in), "-o", str(rel_out)],
                         cwd=SMDR_DIR, capture_output=True, text=True, timeout=1800)
    txt = res.stdout
    if res.returncode != 0 or "Higgs self-coupling" not in txt:
        print(f"  SMDR calc_fit failed for {tag}: {res.stderr[-300:]}")
        return None
    vals = {}
    for key, pat in {"g3": r"g3 = ([\d.]+)", "g2": r"\bg = ([\d.]+)", "gY": r"gp = ([\d.]+)",
                     "yt": r"yt = ([\d.]+)", "yb": r"yb = ([\d.]+)", "ytau": r"ytau = ([\d.]+)",
                     "lam": r"lambda = ([\d.]+)", "Q": r"Q = ([\d.]+)",
                     "MW_pole": r"MW = ([\d.]+);\s+\(\* pole"}.items():
        m = re.search(pat, txt)
        vals[key] = float(m.group(1)) if m else np.nan
    vals["model_file"] = str(out_model)
    cache.write_text(json.dumps(vals, indent=1))
    return vals


def smdr_lambda_curve(tag: str) -> dict | None:
    """SMDR's own running (full state-of-the-art beta functions) from its fitted model."""
    work = cfg.CACHE_DIR / "smdr_runs"
    model = work / f"{tag}_model.dat"
    cache = work / f"{tag}_RGrun.npz"
    if cache.exists():
        return dict(np.load(cache))
    if not (smdr_available() and model.exists()):
        return None
    outfile = work / f"{tag}_RGrun_vs_Q.dat"
    subprocess.run(["./fig_RGrun_vs_Q", "-i", str(Path("../../smdr_runs") / model.name),
                    "-o", str(Path("../../smdr_runs") / outfile.name)],
                   cwd=SMDR_DIR, capture_output=True, text=True, timeout=1800)
    if not outfile.exists():
        return None
    a = np.loadtxt(outfile)
    res = {"Q": a[:, 0], "g3": a[:, 2], "g2": a[:, 3], "gY": a[:, 4], "yt": a[:, 11],
           "lam": a[:, 20]}
    np.savez(cache, **res)
    return res


# ---------------------------------------------------------------------------
# Everything the report needs, cached
# ---------------------------------------------------------------------------

SUMMARY_CACHE = cfg.CACHE_DIR / "stability.json"
CURVES_CACHE = cfg.CACHE_DIR / "stability_curves.npz"
GRID_CACHE = cfg.CACHE_DIR / "stability_grids.npz"

#: Published benchmarks (Buttazzo et al. 2013 eqs. 4.1-4.3; Andreassen et al. 2018 sec. 9)
PUBLISHED = {
    "buttazzo2013": {"lambda_MPl": -0.0143, "yt_MPl": 0.3825, "g1_MPl": 0.6154, "g2_MPl": 0.5055,
                     "g3_MPl": 0.4873, "log10_Lambda_V": 9.5, "Lambda_lambda_over_V": 2.0,
                     "Lambda_I_over_V": 13.0, "stability_bound_mt_eq_4_2": 171.115,
                     "stability_bound_mt_quoted": 171.53},
    "afs2018": {"mu_star": 3.11e17, "lambda_star": -0.0138, "yt_star": 0.402, "g_star": 0.515,
                "gp_star": 0.460, "log10_rate_GeV4": -773, "log10_tau_years": 161,
                "log10_p_light_cone": -606},
}


def _clean(d):
    return {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in d.items()}


def summary(refresh: bool = False) -> dict:
    if SUMMARY_CACHE.exists() and CURVES_CACHE.exists() and not refresh:
        return json.loads(SUMMARY_CACHE.read_text())
    out = {"published": PUBLISHED, "light_cone_H0^4_volume": light_cone_volume()}
    curves = {}
    # 1. benchmarks with the papers' own inputs and matching
    bm = {}
    for ed in ("buttazzo2013", "afs2018"):
        inp = inputs(ed)
        tr = trajectory(inp)
        bm[ed] = _clean(landmarks(tr) | decay(tr)) | {
            "stability_bound_mt": stability_bound_mt(inp),
            "couplings_at_mt": _clean(rge.match(inp))}
        if ed == "afs2018":
            mu, X = tr["mu"], tr["X"]
            i = int(np.argmin(np.where(mu < M_PL, X[:, 6], 9)))
            bm[ed]["at_mu_star"] = {"yt": float(np.sqrt(X[i, 3])), "g": float(np.sqrt(X[i, 1])),
                                    "gp": float(np.sqrt(0.6 * X[i, 0]))}
    out["benchmarks"] = bm

    # 2. our inputs, both matchings
    off = smdr_offsets("pdg2026")
    out["smdr_offsets"] = off
    out["smdr_fit_pdg2026"] = smdr_fit(inputs("pdg2026"), "pdg2026")
    res = {}
    for ed in ("pdg2026", "pdg2024"):
        inp = inputs(ed)
        for name, o in (("buttazzo", None), ("smdr", off)):
            if name == "smdr" and o is None:
                continue
            tr = trajectory(inp, offsets=o)
            res[f"{ed}/{name}"] = _clean(landmarks(tr) | decay(tr)) | {
                "stability_bound_mt": stability_bound_mt(inp, offsets=o),
                "stability_bound_mt_lambda": stability_bound_mt(inp, criterion="lambda",
                                                                offsets=o),
                "couplings_at_mt": _clean(match(inp, offsets=o))}
            if ed == "pdg2026":
                curves[f"{name}_mu"] = tr["mu"]
                curves[f"{name}_X"] = tr["X"]
    out["results"] = res
    best = "pdg2026/smdr" if off is not None else "pdg2026/buttazzo"
    out["headline"] = best

    # 3. loop orders and matching orders (PDG 2026, headline matching)
    inp = inputs("pdg2026")
    loops = {}
    for n in (1, 2, 3):
        tr = trajectory(inp, loops=n, offsets=off)
        loops[str(n)] = _clean(landmarks(tr))
        curves[f"loops{n}_mu"], curves[f"loops{n}_lam"] = tr["mu"], tr["X"][:, 6]
    out["loop_orders"] = loops
    orders = {}
    for o_ in ("LO", "NLO", "NNLO"):
        tr = trajectory(inp, order=o_)
        orders[o_] = _clean(landmarks(tr)) | {"couplings_at_mt": _clean(rge.match(inp, o_))}
        curves[f"order{o_}_lam"] = tr["X"][:, 6]
    out["matching_orders"] = orders

    # 4. uncertainties: one input at a time, +-1 sigma (and +-3 sigma in M_t for the band)
    sig = {"Mt": cfg.M_TOP[1], "alpha_s": cfg.ALPHA_S_MZ[1], "Mh": cfg.M_HIGGS[1]}
    var = {}
    for k, s in sig.items():
        for n in (-3, -1, 1, 3):
            if k != "Mt" and abs(n) == 3:
                continue
            v = getattr(inp, k) + n * s
            tr = trajectory(_with(inp, **{k: v}), offsets=off)
            var[f"{k}{n:+d}"] = _clean(landmarks(tr) | decay(tr))
            curves[f"var_{k}{n:+d}_lam"] = tr["X"][:, 6]
            curves[f"var_{k}{n:+d}_mu"] = tr["mu"]
    out["variations"] = var
    bounds = {}
    for k in ("alpha_s", "Mh"):
        for n in (-1, 1):
            bounds[f"{k}{n:+d}"] = stability_bound_mt(_with(inp, **{k: getattr(inp, k) + n * sig[k]}),
                                                     offsets=off)
    out["bound_variations"] = bounds
    b0 = res[best]["stability_bound_mt"]
    d_as = 0.5 * abs(bounds["alpha_s+1"] - bounds["alpha_s-1"])
    d_mh = 0.5 * abs(bounds["Mh+1"] - bounds["Mh-1"])
    d_match = abs(res["pdg2026/smdr"]["stability_bound_mt"] - res["pdg2026/buttazzo"]
                  ["stability_bound_mt"]) if off is not None else 0.0
    d_higher = 0.15                       # Buttazzo et al.'s higher-order estimate
    d_mc = 0.5                            # Monte Carlo vs pole top mass (indicative)
    gap = inp.Mt - b0
    s_exp = np.sqrt(cfg.M_TOP[1] ** 2 + d_as**2 + d_mh**2)
    s_th = np.sqrt(d_match**2 + d_higher**2)
    out["stability_gap"] = {"Mt": inp.Mt, "bound": b0, "gap_GeV": gap,
                            "sigma_mt": cfg.M_TOP[1], "sigma_alpha_s": d_as, "sigma_mh": d_mh,
                            "sigma_matching": d_match, "sigma_higher_orders": d_higher,
                            "sigma_mc_vs_pole": d_mc,
                            "n_sigma_exp": gap / s_exp,
                            "n_sigma_exp_th": gap / np.hypot(s_exp, s_th),
                            "n_sigma_with_mc_mass": gap / np.sqrt(s_exp**2 + s_th**2 + d_mc**2)}
    smdr_curve = smdr_lambda_curve("pdg2026")
    if smdr_curve is not None:
        curves["smdr_Q"], curves["smdr_lam"] = smdr_curve["Q"], smdr_curve["lam"]
    SUMMARY_CACHE.write_text(json.dumps(out, indent=1, default=float))
    np.savez(CURVES_CACHE, **curves)
    return out


def grids(refresh: bool = False) -> dict:
    """Phase diagrams: a wide view and a zoom on the measured masses."""
    if GRID_CACHE.exists() and not refresh:
        return dict(np.load(GRID_CACHE))
    off = smdr_offsets("pdg2026")
    base = inputs("pdg2026")
    out = {}
    with np.errstate(all="ignore"):
        wide = phase_grid(np.arange(60.0, 220.01, 2.0), np.arange(50.0, 200.01, 2.0), base,
                          offsets=off)
        zoom = phase_grid(np.arange(166.0, 180.01, 0.2), np.arange(118.0, 132.01, 0.2), base,
                          offsets=off)
    for tag, g in (("wide", wide), ("zoom", zoom)):
        for k, v in g.items():
            out[f"{tag}_{k}"] = v
    np.savez(GRID_CACHE, **out)
    return out
