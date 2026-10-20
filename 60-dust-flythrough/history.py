"""The Local Bubble's history: Zucker et al. (2022)'s expanding-shell fit, redone from their data.

Model (their eq. 1, from El-Badry et al. 2019's simulations of clustered
supernovae): the shell radius t Myr after the first explosion is

    R(t) = 83 pc (1 - theta)^(1/5) (E_SN / 1e51 erg)^(1/5) (dt_SNe / 0.1 Myr)^(-1/5)
           (n0 / 1 cm^-3)^(-1/5) (t / 1 Myr)^(3/5),

with cooling efficiency theta = 0.7 and E_SN = 1e51 erg fixed. Its centre is
the mean position of the Upper Centaurus-Lupus and Lower Centaurus-Crux
clusters (whose supernovae made the bubble) at the time of the first
explosion t_exp, in a frame moving with the Local Standard of Rest. Each
younger cluster in their Extended Data Table 2 is assumed born on the shell:
its distance from the centre at its birth (traced back by Zucker et al. with
galpy) should equal R(t_birth - t_exp) within a Gaussian thickness dR.

Free: t_exp, n0, dt_SNe, dR, with Zucker et al.'s priors. Their result:
t_exp = -14.39 (+0.78 -0.74) Myr, n0 = 2.71 (+1.57 -1.02) cm^-3,
dt_SNe = 1.06 (+0.63 -0.39) Myr, dR = 23.3 (+2.5 -2.3) pc; today R = 165 +- 6
pc and v = 6.7 (+0.5 -0.4) km/s, centre (39, 7, -18) pc.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import norm, truncnorm

import catalogs
import config as C

THETA = 0.7
PC_PER_MYR_KMS = 0.9778          # 1 pc/Myr in km/s


def clusters() -> pd.DataFrame:
    """Extended Data Table 1 (ages) joined to the traceback files."""
    t1 = pd.read_csv(C.CATALOG_DIR / "zucker2022" / "Extended_Data_Table_1.dat", sep="\t")
    tracks = catalogs.zucker2022_tracebacks()

    def key(row):
        return f"{row.cluster}_{row.subgroup}".replace(" - ", "_").replace(",", "").replace("/", "_") \
            .replace(" ", "_").replace("-", "_")

    keys = {k.lower().replace("-", "_"): k for k in tracks}
    names = []
    for row in t1.itertuples():
        k = key(row).lower()
        k = {"sco_cen_lower_centaurus_crux_lcc": "sco_cen_lower_centaurus_crux_lcc",
             "sco_cen_upper_centaurus_lupus_ucl": "sco_cen_upper_centaurus_lupus_ucl",
             "sco_cen_upper_scorpius_usco": "sco_cen_upper_scorpius_usco",
             "ophiuchus_rho_oph_population_i": "ophiuchus_rho_oph_population_i",
             "ophiuchus_rho_oph_population_ii": "ophiuchus_rho_oph_population_ii",
             "orion_ngc_2023_2024": "orion_ngc_2023_2024",
             "orion_ngc2068_2071": "orion_ngc2068_2071"}.get(k, k)
        match = [v for kk, v in keys.items() if kk == k or kk.replace("_", "") == k.replace("_", "")]
        names.append(match[0] if match else None)
    t1["track"] = names
    return t1


def position_at(track: pd.DataFrame, t: float) -> np.ndarray:
    return np.array([np.interp(t, track["t"], track[c]) for c in ("x", "y", "z")])


def radius(t_myr, n0, dt_sne):
    t_myr = np.asarray(t_myr, dtype=float)
    return 83.0 * (1 - THETA) ** 0.2 * (dt_sne / 0.1) ** -0.2 * n0 ** -0.2 * np.clip(t_myr, 0, None) ** 0.6


class ExpansionModel:
    names = ["t_exp", "n0", "dt_sne", "dR"]

    def __init__(self):
        cl = clusters()
        tracks = catalogs.zucker2022_tracebacks()
        self.ucl = tracks[cl.loc[cl["subgroup"].str.contains("UCL"), "track"].iloc[0]]
        self.lcc = tracks[cl.loc[cl["subgroup"].str.contains("LCC"), "track"].iloc[0]]
        # the clusters born on the shell: Extended Data Table 2's rows except UCL and LCC
        t2 = pd.read_csv(C.CATALOG_DIR / "zucker2022" / "Extended_Data_Table_2.dat", sep="\t")
        born = []
        for row in t2.itertuples():
            name = row.cluster
            if "UCL" in name or "LCC" in name:
                continue
            region, sub = name.split("(", 1)
            sub = sub.rstrip(")")
            hit = cl[(cl["cluster"].str.strip() == region.strip()) & (cl["subgroup"].str.strip() == sub.strip())]
            if hit.empty:
                raise ValueError(f"no Table 1 row for {name}")
            r = hit.iloc[0]
            born.append({"name": name, "age": float(r["age"]),
                         "birth": position_at(tracks[r["track"]], -float(r["age"]))})
        self.born = born
        self.t_birth = np.array([-b["age"] for b in born])
        self.p_birth = np.array([b["birth"] for b in born])

    def centre(self, t_exp: float) -> np.ndarray:
        return 0.5 * (position_at(self.ucl, t_exp) + position_at(self.lcc, t_exp))

    def loglike(self, th) -> float:
        t_exp, n0, dt_sne, dR = th
        t_i = self.t_birth - t_exp                     # Myr after the first explosion
        if np.any(t_i <= 0):
            return -1e300
        r_i = np.linalg.norm(self.p_birth - self.centre(t_exp), axis=1)
        R_i = radius(t_i, n0, dt_sne)
        return float(-0.5 * np.sum(np.log(2 * np.pi * dR**2) + (R_i - r_i) ** 2 / dR**2))

    @staticmethod
    def prior_transform(u):
        """Zucker et al.: t_exp ~ N(-13, 1) on [-16, -10]; n0 ~ logN(2, x2) on [0.1, 10];
        dt_SNe ~ logN(0.8 Myr, x2) on [0.05, 3]; dR ~ N(15, 5) on [0, 30] pc."""
        def tn(u, mu, sd, lo, hi):
            return truncnorm.ppf(u, (lo - mu) / sd, (hi - mu) / sd, loc=mu, scale=sd)
        ln2 = np.log(2)
        return np.array([tn(u[0], -13, 1, -16, -10),
                         np.exp(tn(u[1], np.log(2), ln2, np.log(0.1), np.log(10))),
                         np.exp(tn(u[2], np.log(0.8), ln2, np.log(0.05), np.log(3))),
                         tn(u[3], 15, 5, 0, 30)])


def fit(force: bool = False, seed: int = 20221) -> dict:
    path = C.CACHE_DIR / "bubble_history.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    import dynesty
    from dynesty.utils import resample_equal
    m = ExpansionModel()
    rng = np.random.default_rng(seed)
    sampler = dynesty.NestedSampler(m.loglike, m.prior_transform, 4, nlive=800, rstate=rng)
    sampler.run_nested(dlogz=0.1, print_progress=False)
    res = sampler.results
    s = resample_equal(res.samples, res.importance_weights(), rstate=rng)
    t_now = -s[:, 0]
    R_now = radius(t_now, s[:, 1], s[:, 2])
    v_now = 0.6 * R_now / t_now * PC_PER_MYR_KMS
    centres = np.array([m.centre(t) for t in s[:400, 0]])
    pct = lambda a: [float(v) for v in np.percentile(a, [16, 50, 84])]
    out = {"n_clusters": len(m.born), "clusters": [b["name"] for b in m.born],
           "t_exp_myr": pct(s[:, 0]), "n0_cm3": pct(s[:, 1]), "dt_sne_myr": pct(s[:, 2]), "dR_pc": pct(s[:, 3]),
           "R_now_pc": pct(R_now), "v_now_kms": pct(v_now), "n_sne": pct(t_now / s[:, 2]),
           "centre_pc": [float(v) for v in np.median(centres, axis=0)], "logz": float(res.logz[-1])}
    path.write_text(json.dumps(out, indent=1))
    np.save(C.CACHE_DIR / "bubble_history_samples.npy", s)
    return out


def shell_vs_model(model_centre, R_model) -> dict:
    """Distance of our dust wall (each sightline's first peak) from the model's centre."""
    import bubble
    shell = bubble.shell_model()
    ok = shell["d_peak"].notna() & (shell["d_peak"] < 400)
    P = shell.loc[ok, ["x", "y", "z"]].to_numpy()
    r = np.linalg.norm(P - np.asarray(model_centre), axis=1)
    return {"median_r_pc": float(np.median(r)), "p16_p84_pc": [float(np.quantile(r, 0.16)), float(np.quantile(r, 0.84))],
            "frac_within_model_pm_30pc": float(np.mean(np.abs(r - R_model) < 30)), "R_model_pc": float(R_model)}


if __name__ == "__main__":
    res = fit(force=True)
    print(json.dumps(res, indent=1))
    print(json.dumps(shell_vs_model(res["centre_pc"], res["R_now_pc"][1]), indent=1))
