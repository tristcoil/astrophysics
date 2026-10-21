"""Is GW170817-like europium production enough for the whole Milky Way?

1. Solar r-process matter, element by element: Asplund et al. (2009)'s solar
   abundances times (1 - s fraction) from Bisterzo et al. (2014). This gives
   the share of europium in r-process material, for the whole range Bisterzo
   covers (krypton upward) and for the lanthanides and heavier.
2. The kilonova's europium: ejecta mass x europium share, for three readings
   of which components carry heavy r-process material.
3. The Galaxy's europium: (stars + gas) x the solar europium mass fraction x
   the mean europium level of GALAH's disc stars, times Eu's r-process share.
4. Mergers needed over the Galaxy's 12 Gyr, against mergers available: the
   GWTC-4.0 local rate scaled to the Milky Way by stellar mass (with project
   J3's local stellar mass density), and to the Galaxy's history by the
   cosmic star-formation history of Madau & Dickinson (2014) convolved with a
   1/t delay-time distribution starting at 20 Myr.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import (J3_LOG_RHO_STAR, LIT_DIR, MW_AGE_GYR, MW_GAS_MASS, MW_STELLAR_MASS,
                    REFERENCE)

#: Standard atomic weights (IUPAC), for the elements used here.
ATOMIC_WEIGHT = {
    "H": 1.008, "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Ru": 101.07, "Rh": 102.906, "Pd": 106.42, "Ag": 107.868,
    "Cd": 112.414, "In": 114.818, "Sn": 118.710, "Sb": 121.760, "Te": 127.60, "I": 126.904,
    "Xe": 131.293, "Cs": 132.905, "Ba": 137.327, "La": 138.905, "Ce": 140.116, "Pr": 140.908,
    "Nd": 144.242, "Sm": 150.36, "Eu": 151.964, "Gd": 157.25, "Tb": 158.925, "Dy": 162.500,
    "Ho": 164.930, "Er": 167.259, "Tm": 168.934, "Yb": 173.045, "Lu": 174.967, "Hf": 178.49,
    "Ta": 180.948, "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.967, "Hg": 200.592, "Tl": 204.38, "Pb": 207.2, "Bi": 208.980, "Th": 232.038,
    "U": 238.029,
}
#: Asplund et al. (2009), Table 4: present-day photospheric hydrogen mass fraction.
X_H_SUN = 0.7381
NOBLE = {"Kr", "Xe"}          # depleted in meteorites: use the (indirect) solar values
LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
               "Yb", "Lu"]


def solar_table() -> pd.DataFrame:
    t = pd.read_csv(LIT_DIR / "asplund2009_table1.csv")
    t["element"] = t["element"].str.strip()
    use_met = ~t["element"].isin(NOBLE) & t["meteorites"].notna() & (t["Z"] >= 31)
    t["logeps"] = np.where(use_met, t["meteorites"], t["photosphere"])
    return t.set_index("element")


def s_fractions() -> pd.Series:
    b = pd.read_csv(LIT_DIR / "bisterzo2014_sprocess.csv")
    iso = b["Isotope"].astype(str).str.strip()
    elem = b[iso.str.fullmatch(r"[A-Z][a-z]?")]
    return pd.Series(elem["TW+"].to_numpy(float) / 100.0, index=iso[elem.index])


def solar_r_process() -> pd.DataFrame:
    """Mass fraction of each element in the Sun, and of its r-process part
    (1 - s; thorium and uranium are pure r)."""
    sol = solar_table()
    s = s_fractions()
    rows = []
    for el, w in ATOMIC_WEIGHT.items():
        if el == "H" or el not in sol.index or not np.isfinite(sol.loc[el, "logeps"]):
            continue
        x = X_H_SUN * 10 ** (sol.loc[el, "logeps"] - 12.0) * w / ATOMIC_WEIGHT["H"]
        s_frac = 0.0 if el in ("Th", "U") else s.get(el, np.nan)
        if not np.isfinite(s_frac):
            continue
        rows.append(dict(element=el, Z=int(sol.loc[el, "Z"]), x_sun=x, s_frac=float(s_frac),
                         x_r=x * (1 - min(s_frac, 1.0))))
    return pd.DataFrame(rows).set_index("element")


def europium_shares() -> dict:
    r = solar_r_process()
    eu = r.loc["Eu", "x_r"]
    heavy = r[r["Z"] >= 57]["x_r"].sum()
    lan = r[r.index.isin(LANTHANIDES)]["x_r"].sum()
    return dict(x_eu_sun=float(r.loc["Eu", "x_sun"]), eu_r_fraction=float(1 - r.loc["Eu", "s_frac"]),
                x_r_kr_up=float(r["x_r"].sum()), x_r_la_up=float(heavy),
                eu_share_kr_up=float(eu / r["x_r"].sum()), eu_share_la_up=float(eu / heavy),
                lanthanide_share_kr_up=float(lan / r["x_r"].sum()),
                n_elements=int(len(r)))


def kilonova_europium(m_blue, m_purple, m_red, shares) -> dict:
    """Europium made by the kilonova, for three readings of the ejecta:

    * 'all': every component carries the solar r-process pattern from Kr up;
    * 'purple+red': only the lanthanide-bearing (kappa >= 3) components do;
      the fast blue ejecta made only light r-process elements, no Eu;
    * 'red only': only the lanthanide-rich (kappa = 10) component does.
    """
    f = shares["eu_share_kr_up"]
    return {"all": (m_blue + m_purple + m_red) * f, "purple+red": (m_purple + m_red) * f,
            "red only": m_red * f}


def galaxy_europium(mean_eu_h_linear: float) -> dict:
    """r-process europium in the Milky Way's stars and gas today (Msun)."""
    sh = europium_shares()
    m_star, m_star_err = MW_STELLAR_MASS
    m = (m_star + MW_GAS_MASS) * sh["x_eu_sun"] * mean_eu_h_linear * sh["eu_r_fraction"]
    frac_err = m_star_err / (m_star + MW_GAS_MASS)
    return dict(m_eu=float(m), m_eu_err=float(m * frac_err), x_eu_sun=sh["x_eu_sun"])


# ---------------------------------------------------------------------------
# Merger rates through cosmic time
# ---------------------------------------------------------------------------


def madau_dickinson(z):
    """Cosmic star-formation rate density (Msun/yr/Mpc^3), Madau & Dickinson 2014."""
    return 0.015 * (1 + z) ** 2.7 / (1 + ((1 + z) / 2.9) ** 5.6)


def history_factor(*, t_min_gyr=0.02, span_gyr=MW_AGE_GYR, slope=-1.0):
    """Mean merger rate over the last ``span_gyr`` divided by today's, for
    mergers that follow star formation with delays ~ t^slope (t > t_min)."""
    from astropy.cosmology import FlatLambdaCDM

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)     # Madau & Dickinson's cosmology
    z = np.logspace(-4, np.log10(20), 4000)[::-1]
    t = cosmo.age(z).to("Gyr").value           # increasing
    psi = madau_dickinson(z)
    grid = np.linspace(t[0], t[-1], 6000)
    psi_g = np.interp(grid, t, psi)
    dt = grid[1] - grid[0]
    tau = np.arange(1, grid.size) * dt
    dtd = np.where(tau > t_min_gyr, tau ** slope, 0.0)
    rate = np.array([np.sum(psi_g[:i][::-1] * dtd[:i]) * dt for i in range(1, grid.size)])
    tt = grid[1:]
    now = rate[-1]
    window = tt > tt[-1] - span_gyr
    return float(np.mean(rate[window]) / now), tt, rate / now


def merger_budget(m_eu_per_event: float, m_eu_galaxy: float, r0_gpc3_yr) -> dict:
    """Mergers the Milky Way needed versus mergers it had."""
    n_needed = m_eu_galaxy / m_eu_per_event
    rate_needed = n_needed / (MW_AGE_GYR * 1e3)          # per Myr, averaged
    per_mass = 1e-9 / 10 ** J3_LOG_RHO_STAR             # mergers / yr / Msun for R0 = 1
    rate_now = np.asarray(r0_gpc3_yr) * per_mass * MW_STELLAR_MASS[0] * 1e6   # per Myr
    f_hist, _, _ = history_factor()
    return dict(n_needed=float(n_needed), rate_needed_per_myr=float(rate_needed),
                rate_now_per_myr=rate_now.tolist(), history_factor=f_hist,
                rate_avg_per_myr=(rate_now * f_hist).tolist(),
                supply_over_demand=(rate_now * f_hist / rate_needed).tolist())


def villar_2017_check(shares) -> dict:
    """Villar et al. (2017) estimated the r-process mass each merger must make
    from Milky Way production rates (Qian 2000; Bauswein et al. 2014) and the
    2017 rate; repeat their arithmetic with today's rate."""
    v = REFERENCE["villar2017"]
    scale = REFERENCE["bns_rate_2017"] / np.sqrt(np.prod(REFERENCE["gwtc4_bns_rate"]))
    return dict(needed_light_2017=v["needed_light_msun"], needed_heavy_2017=v["needed_heavy_msun"],
                needed_light_now=v["needed_light_msun"] * scale,
                needed_heavy_now=v["needed_heavy_msun"] * scale, rate_scale=float(scale))
