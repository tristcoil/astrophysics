"""Removing the solar corona from Cassini's Doppler with two radio frequencies.

Free electrons advance a radio wave's phase by 40.3 N_e / f^2 per metre, so
the corona adds to each leg of the link a fractional frequency shift
A/f^2, where A is proportional to the rate of change of the electron
column along that leg. For a two-way link with uplink frequency f_u and
downlink f_d (the transponder preserves fractional shifts),

    y_link = y_0 + A_u / f_u^2 + A_d / f_d^2,

with y_0 everything non-dispersive: gravity, the Shapiro rate, trajectory
errors, the troposphere. The uplink (sent at t1) and the downlink (received
at t3) cross the corona about two hours apart, so A_u and A_d differ.

**Three links** (Bertotti, Comoretto & Iess 1993, A&A 269, 608): Cassini
could carry X/X (7.2 up, 8.4 GHz down), X/Ka (7.2 up, 32.0 down) and Ka/Ka
(34.3 up, 32.0 down) at once. X/X and X/Ka share the uplink, X/Ka and Ka/Ka
(nearly) the downlink frequency, so three equations give y_0, A_u and A_d
exactly (``three_link``): y_0 = 1.0457 y_KaKa + 0.0287 y_XKa - 0.0744 y_XX.
The X/Ka link, however, was recorded only by the open-loop Radio Science
Receivers (RSR files, ~120 MB a pass): the closed-loop data in the archive
(ODFs and TDFs) do not contain it, because DSS 25's single Ka-band
closed-loop channel tracked the Ka/Ka link (apart from ~1 h on June 8 when it
held the X/Ka downlink, labelled Ka/Ka and screened out; ``cassini_tdf``).

**Two links** (what the closed-loop data allow): X/X and Ka/Ka at DSS 25,
simultaneous on 23 of the 29 days with two-way data at DSS 25 (not June 18-20,
when the Ka uplink was off). With P = 1/f_u^2 + 1/f_d^2 for each link,

    y_0 = (P_X y_Ka - P_Ka y_X) / (P_X - P_Ka) = 1.0576 y_Ka - 0.0576 y_X,

exact when A_u = A_d. That holds for the slowly changing corona (the mean
density gradient the line of sight sweeps through), so the systematic
plasma signal, the part that could mimic or mask the Shapiro rate, cancels.
For turbulence that differs between the legs the residue is
2.2e-4 (A_d - A_u) in GHz^-2 units: about 4 times less than in Ka/Ka alone
(1.3e-3 rms) and 77 times less than in X/X.

The scale of the cancellation, for the record: P_X = 0.0335 GHz^-2,
P_Ka = 0.00182 GHz^-2, so X/X carries 18.4 times Ka/Ka's plasma.

Found in the data: away from the Sun (b > 8 R_sun) the combination is
quieter than Ka/Ka alone (60-s noise 4e-14 against 6e-14 at b = 12-20), but
inside 5 R_sun it is noisier (1.3e-12 against 7.6e-13 at b = 3-5): there the
closed-loop X/X data carry noise that does not scale as 1/f^2 (receiver
tracking of a strongly scintillating carrier), which the -0.0576 weight
passes on. The closest days therefore add little to the fit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import cassini_doppler as D


def leg_frequencies(ul: int, dl: int, f_t) -> tuple[np.ndarray, np.ndarray]:
    """Uplink and downlink sky frequencies (Hz) of a link."""
    f_t = np.asarray(f_t, dtype=float)
    return f_t, D.M2[(ul, dl)] * f_t


def plasma_weight(ul: int, dl: int, f_t) -> np.ndarray:
    """P = 1/f_u^2 + 1/f_d^2 in GHz^-2."""
    fu, fd = leg_frequencies(ul, dl, f_t)
    return (1e9 / fu) ** 2 + (1e9 / fd) ** 2


#: label and pass-number offset of the combined link
COMBO = "X/X+Ka/Ka"
COMBO_PASS_OFFSET = 100000


def two_link_1s(cl: pd.DataFrame, station: int = 25) -> pd.DataFrame:
    """The plasma-free combination at the 1-s level, from cleaned 1-s points
    of X/X and Ka/Ka with the same time tag (so its white noise comes from
    its own scatter: the two links' plasma noise is correlated). Geometry,
    y_GR and partials are the Ka/Ka point's (same legs, same times). Its
    pass and track are the Ka/Ka pass's."""
    x = cl[(cl["station"] == station) & (cl["ul_band"] == 2) & (cl["dl_band"] == 2)]
    k = cl[(cl["station"] == station) & (cl["ul_band"] == 3) & (cl["dl_band"] == 3)]
    m = k.merge(x[["t1950", "y_obs", "f_t"]], on="t1950", suffixes=("", "_x"))
    p_x = plasma_weight(2, 2, m["f_t_x"])
    p_k = plasma_weight(3, 3, m["f_t"])
    w_x, w_k = -p_k / (p_x - p_k), p_x / (p_x - p_k)
    out = m.drop(columns=["y_obs_x", "f_t_x"]).copy()
    out["y_obs"] = w_k * m["y_obs"] + w_x * m["y_obs_x"]
    out["y_plasma_x"] = (m["y_obs_x"] - m["y_obs"]) / (p_x - p_k) * p_x
    out["link"] = COMBO
    out["ul_band"], out["dl_band"] = 0, 0
    out["pass_id"] = out["pass_id"] + COMBO_PASS_OFFSET
    return out


def two_link(np_: pd.DataFrame, station: int = 25) -> pd.DataFrame:
    """Plasma-free combination of simultaneous X/X and Ka/Ka *normal points*
    (a check of ``two_link_1s``; its white error treats the two links' noise
    as independent, which overstates it where plasma dominates both).

    Returns one row per common time bin with the combined residual y_obs,
    GR's signal, a white-noise error, the estimated X-band plasma
    contribution, and the geometry/partials of the Ka/Ka point (identical
    up and down legs, so identical y_GR and partials to 1e-16)."""
    x = np_[(np_["station"] == station) & (np_["link"] == "X/X")]
    k = np_[(np_["station"] == station) & (np_["link"] == "Ka/Ka")]
    m = k.merge(x[["tbin", "y_obs", "y_err_white", "f_t", "n", "y_gr"]], on="tbin", suffixes=("", "_x"))
    p_x = plasma_weight(2, 2, m["f_t_x"])
    p_k = plasma_weight(3, 3, m["f_t"])
    w_x = -p_k / (p_x - p_k)
    w_k = p_x / (p_x - p_k)
    out = m.copy()
    out["y_obs_ka"] = m["y_obs"]
    out["y_obs_x"] = m["y_obs_x"]
    out["y_obs"] = w_k * m["y_obs"] + w_x * m["y_obs_x"]
    out["y_err_white"] = np.sqrt((w_k * m["y_err_white"]) ** 2 + (w_x * m["y_err_white_x"]) ** 2)
    a = (m["y_obs_x"] - m["y_obs"]) / (p_x - p_k)             # common plasma A (GHz^2 units)
    out["y_plasma_x"] = a * p_x
    out["y_plasma_ka"] = a * p_k
    out["w_x"], out["w_ka"] = w_x, w_k
    out["link"] = "X/X+Ka/Ka"
    out["gr_mismatch"] = m["y_gr"] - m["y_gr_x"]
    return out


def three_link(y_xx, y_xka, y_kaka, f_ux, f_uka) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact separation (y_0, A_u at X, A_d) from X/X, X/Ka and Ka/Ka
    fractional residuals (for open-loop data; not used with the archive's
    closed-loop data, which lack X/Ka). Frequencies in Hz."""
    f_ux, f_uka = np.asarray(f_ux, float), np.asarray(f_uka, float)
    a = np.zeros(np.shape(f_ux) + (3, 3))
    a[..., 0, :] = np.stack([np.ones_like(f_ux), (1e9 / f_ux) ** 2, (1e9 / (D.M2[(2, 2)] * f_ux)) ** 2], -1)
    a[..., 1, :] = np.stack([np.ones_like(f_ux), (1e9 / f_ux) ** 2, (1e9 / (D.M2[(2, 3)] * f_ux)) ** 2], -1)
    a[..., 2, :] = np.stack([np.ones_like(f_ux), (1e9 / f_uka) ** 2, (1e9 / (D.M2[(3, 3)] * f_uka)) ** 2], -1)
    # A_u is per link (X uplink for the first two, Ka uplink for the third):
    # the same electron column, so A_u is common when expressed per 1/f^2
    sol = np.linalg.solve(a, np.stack([y_xx, y_xka, y_kaka], -1)[..., None])[..., 0]
    return sol[..., 0], sol[..., 1], sol[..., 2]


def plasma_template(np_: pd.DataFrame, index: float = 1.3) -> np.ndarray:
    """Shape of the mean coronal contribution to a link's y: for an electron
    density falling as r^-(index+1) the column across the corona varies as
    b^-index, so each leg adds d(b^-index)/dt / f^2 (in units of R_sun and
    GHz^-2 per day). Time derivatives come from the smooth b(t) within each
    pass. For the X/X + Ka/Ka combination the template is the combination of
    the two links' templates (what survives of a b-dependent corona whose
    up and down legs differ). The fitted amplitude carries sign and size."""
    out = np.zeros(len(np_))
    for pid, g in np_.groupby("pass_id"):
        if len(g) < 3:
            continue
        t = g["t_tdb"].to_numpy() / 86400.0
        ul, dl = int(g["ul_band"].iloc[0]), int(g["dl_band"].iloc[0])
        legs = {col: np.gradient(g[col].to_numpy() ** (-index), t) for col in ("b_up_rsun", "b_down_rsun")}
        if (ul, dl) == (0, 0):
            # the X/X + Ka/Ka combination: its weights applied to each link's legs
            fk = g["f_t"].to_numpy()
            fx = np.full(len(g), 7.175e9)
            p_x, p_k = plasma_weight(2, 2, fx), plasma_weight(3, 3, fk)
            w_x, w_k = -p_k / (p_x - p_k), p_x / (p_x - p_k)
            s = 0.0
            for (u, d_, f), w in (((2, 2, fx), w_x), ((3, 3, fk), w_k)):
                fu, fd = leg_frequencies(u, d_, f)
                s = s + w * (legs["b_up_rsun"] * (1e9 / fu) ** 2 + legs["b_down_rsun"] * (1e9 / fd) ** 2)
        else:
            fu, fd = leg_frequencies(ul, dl, g["f_t"].to_numpy())
            s = legs["b_up_rsun"] * (1e9 / fu) ** 2 + legs["b_down_rsun"] * (1e9 / fd) ** 2
        out[np_.index.get_indexer(g.index)] = s
    return out
