"""Every figure of the project, drawn from the cached analysis results."""

from __future__ import annotations

import io

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C
import lightcurves as LC
import literature as lit
import progenitors as P
import rates as R
from plotstyle import (AQUA, AXIS, BLUE, GREEN, GRID, INK, INK_SECONDARY, MAGENTA, MUTED, ORANGE,
                       RED, SURFACE, VIOLET, YELLOW, annotate, caption, save, use_style)

use_style()

FAMILY_COLOR = {"Ia": BLUE, "II": ORANGE, "Ibc": AQUA, "IIn": YELLOW, "SLSN": MAGENTA, "TDE": GREEN}
GREY_CLOUD = "#c9c7bf"
SOURCE_BTS = ("Data: ZTF Bright Transient Survey (Fremling et al. 2020; Perley et al. 2020), "
              "catalogue snapshot of 4 Oct 2026; statistical sample = quality + purity cuts, "
              "peak brighter than 18.5 mag, peak between 1 Jun 2018 and 31 Oct 2025.")


def _out(name: str) -> str:
    return str(C.FIGURE_DIR / name)


def _stat(df):
    import bts
    return df[bts.statistical_sample(df) & df["family"].notna() & df["M_peak"].notna()]


# ---------------------------------------------------------------------------
# 1. The zoo in luminosity-timescale space
# ---------------------------------------------------------------------------

def fig_zoo(res, df) -> list:
    s = _stat(df)
    s = s[~s["duration_is_limit"] & s["duration"].notna()].copy()
    s["dur_rest"] = s["duration"] / (1 + s["z_cmb"])
    fams = list(C.FAMILY_ORDER)
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), sharex=True, sharey=True)
    for ax, fam in zip(axes.flat, fams):
        ax.scatter(s["dur_rest"], s["M_peak"], s=3, color=GREY_CLOUD, lw=0, rasterized=True)
        f = s[s["family"] == fam]
        ax.scatter(f["dur_rest"], f["M_peak"], s=7, color=FAMILY_COLOR[fam], lw=0.3,
                   edgecolor=SURFACE, rasterized=True)
        ax.set_title(C.FAMILY_LABEL[fam], fontsize=10)
        ax.text(0.97, 0.04, f"{len(f):,} events", transform=ax.transAxes, ha="right",
                fontsize=8, color=INK_SECONDARY)
        ax.set_xscale("log")
    ax0 = axes[0, 0]
    ax0.set_ylim(-12.5, -23.5)
    ax0.set_xlim(4, 400)
    for ax in axes[1]:
        ax.set_xlabel("time above half of peak brightness (rest-frame days)")
    for ax in axes[:, 0]:
        ax.set_ylabel("peak absolute magnitude")
    # landmarks
    famous = {"SN2023ixf": ("SN 2023ixf (M101)", (2.2, 1.6)),
              "SN2020jfo": ("SN 2020jfo (M61)", (2.0, -1.4)),
              "SN2022hrs": ("SN 2022hrs (NGC 4647)", (1.7, -1.3))}
    for name, (label, (fx, dy)) in famous.items():
        row = df[df["iau_id"] == name]
        if len(row) and np.isfinite(row["M_peak"].iat[0]) and np.isfinite(row["duration"].iat[0]):
            x = row["duration"].iat[0] / (1 + row["z_cmb"].iat[0])
            fam = row["family"].iat[0]
            ax = axes.flat[fams.index(fam)] if fam in fams else axes.flat[0]
            ax.scatter([x], [row["M_peak"].iat[0]], s=40, facecolor="none", edgecolor=INK, lw=0.9)
            annotate(ax, label, (x, row["M_peak"].iat[0]), (x * fx, row["M_peak"].iat[0] + dy))
    fig.suptitle("The supernova zoo: how bright, and for how long", x=0.06, ha="left")
    fig.tight_layout()
    caption(fig, SOURCE_BTS + " Grey: all classified transients; colour: the family named. "
            "Absolute magnitudes in the band of the observed peak, corrected for Milky Way dust but "
            "not for dust in the host galaxy. Durations are the survey's own (rise + fade to half "
            "flux) divided by 1+z; events whose decline was not followed to half flux are left out.")
    return save(fig, _out("01_supernova_zoo"))


# ---------------------------------------------------------------------------
# 2. Light-curve gallery
# ---------------------------------------------------------------------------

GALLERY = [("SN Ia", "Type Ia", "a white dwarf detonates"),
           ("SN Ia-91bg", "Type Ia-91bg", "faint, fast thermonuclear"),
           ("SN Iax", "Type Iax", "a failed / partial Ia"),
           ("TDE", "Tidal disruption", "a star torn apart (not a SN)"),
           ("SN II", "Type II", "red supergiant: the plateau"),
           ("SN IIb", "Type IIb", "a thin hydrogen layer left"),
           ("SN IIn", "Type IIn", "crashing into old ejecta"),
           ("SN Ibn", "Type Ibn", "helium-rich shells, fast"),
           ("SN Ib", "Type Ib", "stripped of hydrogen"),
           ("SN Ic", "Type Ic", "stripped of H and He"),
           ("SN Ic-BL", "Type Ic-BL", "broad lines: 20,000+ km/s"),
           ("SLSN-I", "Superluminous I", "10-100x brighter, slow")]


def gallery_picks(df, per_class=4) -> dict:
    """For each class, the gallery objects with the most r-band nights around peak."""
    out = {}
    have = {p.stem for p in C.ALERCE_DIR.glob("*.json")}
    for tns, _, _ in GALLERY:
        cand = df[(df["tns_type"] == tns) & df["ztf_id"].isin(have) & df["mu"].notna()
                  & df["quality"]]
        scored = []
        for _, row in cand.iterrows():
            got = LC.load_lc(row["ztf_id"])
            if got is None:
                continue
            nl = LC.nightly(got[0])
            ph = (nl["mjd"] + 2400000.5 - row["peak_jd"]) / (1 + row["z_cmb"])
            n = int(((nl["fid"] == 2) & (ph > -20) & (ph < 60)).sum())
            scored.append((n, row["ztf_id"]))
        scored.sort(reverse=True)
        out[tns] = [z for _, z in scored[:per_class]]
    return out


def fig_gallery(res, df) -> list:
    picks = gallery_picks(df)
    fig, axes = plt.subplots(3, 4, figsize=(12, 8.2), sharey=True)
    dfi = df.set_index("ztf_id")
    for ax, (tns, title, blurb) in zip(axes.flat, GALLERY):
        fam = C.FAMILY.get(C.CLASS_MAP.get(tns, ""), "TDE")
        col = FAMILY_COLOR.get(fam, MUTED)
        xmax = 220 if tns in ("SLSN-I", "TDE") else 130
        for k, oid in enumerate(picks.get(tns, [])):
            row = dfi.loc[oid]
            det, _ = LC.load_lc(oid)
            nl = LC.nightly(det)
            z = row["z_cmb"]
            ph = (nl["mjd"] + 2400000.5 - row["peak_jd"]) / (1 + z)
            r = (nl["fid"] == 2) & (ph > -40) & (ph < xmax)
            M = nl["mag"] - row["mu"] - C.A_R_PER_AV * row["A_V"] + 2.5 * np.log10(1 + z)
            alpha = 1.0 if k == 0 else 0.45
            xs, ys = ph[r].to_numpy(), M[r].to_numpy()
            gap = np.concatenate([[False], np.diff(xs) > 25])
            xs_l = np.insert(xs, np.flatnonzero(gap), np.nan)
            ys_l = np.insert(ys, np.flatnonzero(gap), np.nan)
            ax.plot(xs_l, ys_l, "-", color=col, lw=1.0, alpha=alpha)
            ax.plot(ph[r], M[r], "o", color=col, ms=2.4, alpha=alpha, mec="none")
            if k == 0:
                name = row["iau_id"] if row["iau_id"] != "-" else oid
                ax.text(0.97, 0.05, f"{name}  z = {row['redshift']:.3f}", transform=ax.transAxes,
                        ha="right", fontsize=7, color=INK_SECONDARY)
        ax.set_title(title, fontsize=10)
        ax.text(0.03, 0.93, blurb, transform=ax.transAxes, fontsize=7.5, color=INK_SECONDARY,
                va="top")
        ax.set_xlim(-40, xmax)
    axes[0, 0].set_ylim(-13.0, -22.8)
    for ax in axes[2]:
        ax.set_xlabel("rest-frame days from peak")
    for ax in axes[:, 0]:
        ax.set_ylabel("absolute magnitude (ZTF r)")
    fig.suptitle("Twelve ways to explode: ZTF light curves of each class", x=0.06, ha="left")
    fig.tight_layout()
    caption(fig, "Data: ZTF public alert photometry via the ALeRCE broker (Förster et al. 2021); "
            "classifications from TNS via the BTS. In each panel the four events of that class "
            "with the best r-band coverage around peak (the first drawn solid); nightly means. "
            "Alert photometry stops near 20.5 mag, so faint tails are cut off.")
    return save(fig, _out("02_light_curve_gallery"))


# ---------------------------------------------------------------------------
# 3. Sky map
# ---------------------------------------------------------------------------

def fig_sky(res, df) -> list:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    s = _stat(df)
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(111, projection="mollweide")

    def lonlat(ra, dec):
        lon = -np.radians(((np.asarray(ra) + 180) % 360) - 180)   # RA increases to the left
        return lon, np.radians(dec)

    # shade the parts of the sky the statistical sample cannot contain
    ra_g, dec_g = np.meshgrid(np.linspace(0, 360, 721), np.linspace(-90, 90, 361))
    b = SkyCoord(ra=ra_g.ravel() * u.deg, dec=dec_g.ravel() * u.deg).galactic.b.deg.reshape(ra_g.shape)
    excluded = (np.abs(b) < C.FOOTPRINT_MIN_ABS_B) | (dec_g < C.FOOTPRINT_MIN_DEC)
    lon, lat = lonlat(ra_g, dec_g)
    order = np.argsort(lon[0])
    cs = ax.contourf(lon[:, order], lat[:, order], excluded[:, order].astype(float), levels=[0.5, 1.5],
                     colors=["#ebe9e2"], zorder=0)
    cs.set_rasterized(True)   # the vector outline of the band would weigh ~6 MB in the SVG
    groups = [("other", s[~s["family"].isin(["Ia"]) & ~s["is_cc"]], GREY_CLOUD, "other (TDEs, Ia-CSM)"),
              ("Ia", s[s["family"] == "Ia"], BLUE, "Type Ia"),
              ("CC", s[s["is_cc"]], ORANGE, "core collapse")]
    for key, g, col, lab in groups:
        x, y = lonlat(g["ra"], g["dec"])
        ax.scatter(x, y, s=3.5 if key != "other" else 5, color=col, lw=0, label=f"{lab} ({len(g):,})",
                   rasterized=True)
    ax.set_xticks(np.radians(np.arange(-150, 181, 30)))
    ax.set_xticklabels([f"{(-v) % 360:d}°" for v in range(-150, 181, 30)], fontsize=7, color=MUTED)
    ax.grid(True, color=GRID, lw=0.5)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.17), ncol=3, markerscale=3)
    ax.text(np.radians(-100), np.radians(-50), "below dec −26°:\nnot visible long enough\nfrom Palomar",
            fontsize=7.5, color=INK_SECONDARY, ha="center")
    gp = SkyCoord(l=75 * u.deg, b=0 * u.deg, frame="galactic").icrs
    gx, gy = lonlat(gp.ra.deg, gp.dec.deg)
    ax.text(gx, gy, "Milky Way plane (|b| < 7°)", fontsize=7.5, color=INK_SECONDARY,
            ha="center", va="center", rotation=-62)
    ax.set_title("Seven years of exploding stars over the northern sky", pad=18)
    caption(fig, SOURCE_BTS + " Equatorial coordinates (Mollweide). The shaded band and cap are "
            "excluded by the survey's own cuts. Uneven density elsewhere reflects how long each "
            "part of the sky was observable and how dusty it is.", y=0.02)
    return save(fig, _out("03_sky_map"))


# ---------------------------------------------------------------------------
# 4. Magnitude-limited vs volume-limited
# ---------------------------------------------------------------------------

def fig_fractions(res, df) -> list:
    ml = res["rates"]["mag_limited"]
    ml = ml[ml["m_lim"] >= 15.0]
    fig = plt.figure(figsize=(12, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])
    ax = fig.add_subplot(gs[0])
    x = ml["m_lim"].to_numpy()
    base = np.zeros_like(x)
    for fam in C.FAMILY_ORDER:
        y = ml[fam].to_numpy()
        ax.fill_between(x, base, base + y, color=FAMILY_COLOR[fam], lw=0, alpha=0.9)
        # direct label at the right end, inside the band if it is wide enough
        mid = base[-1] + y[-1] / 2
        if y[-1] > 0.03:
            ax.text(x[-1] + 0.05, mid, f"{C.FAMILY_LABEL[fam].split(' (')[0]}  {100 * y[-1]:.0f}%",
                    va="center", fontsize=8, color=INK)
        base = base + y
    thin = ml["SLSN"].iat[-1] + ml["TDE"].iat[-1]
    ax.text(x[-1] + 0.05, 1.0 - thin / 2 + 0.012, f"superluminous + TDE  {100 * thin:.0f}%",
            va="center", fontsize=8, color=INK)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel("survey depth: peak brighter than ... (apparent mag)")
    ax.set_ylabel("share of classified supernovae")
    ax.set_title("A magnitude-limited survey sees mostly Type Ia...")
    ax.axvline(C.M_LIM, color=INK_SECONDARY, lw=0.8)
    ax.text(C.M_LIM - 0.05, 0.02, "BTS classification\nlimit", ha="right", fontsize=7.5,
            color=INK_SECONDARY)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))

    ax2 = fig.add_subplot(gs[1])
    fr = res["rates"]["main"]["fractions"]["M<-14.0"]
    ml185 = ml.iloc[(ml["m_lim"] - C.M_LIM).abs().argmin()]
    ia_ml = ml185["Ia"]
    ibc_ml = ml185["Ibc"]
    ii_ml = ml185["II"] + ml185["IIn"]
    tot = ia_ml + ibc_ml + ii_ml
    rows = [("BTS: what it sees\n(m < 18.5, this work)", ia_ml / tot, ibc_ml / tot, ii_ml / tot),
            ("LOSS: magnitude-limited\n(Li et al. 2011)", *(lit.LI_2011_MAGLIM[k][0] / 100 for k in ("Ia", "Ibc", "II"))),
            ("BTS: per unit volume\n(1/Vmax, this work)", fr["Ia"]["frac"], fr["Ibc"]["frac"], fr["II"]["frac"]),
            ("LOSS: volume-limited\n(Li et al. 2011)", *(lit.LI_2011_VOLUME[k][0] / 100 for k in ("Ia", "Ibc", "II")))]
    for i, (lab, a, b_, c) in enumerate(rows[::-1]):
        left = 0
        for val, col, name in ((a, BLUE, "Ia"), (b_, AQUA, "Ib/c"), (c, ORANGE, "II")):
            ax2.barh(i, val - 0.004, left=left + 0.002, height=0.62, color=col)
            if val > 0.07:
                ax2.text(left + val / 2, i, f"{name}\n{100 * val:.0f}%", ha="center", va="center",
                         fontsize=7.5, color=SURFACE if col == BLUE else INK)
            left += val
    ax2.set_yticks(range(len(rows)))
    ax2.set_yticklabels([r[0] for r in rows[::-1]], fontsize=8)
    ax2.set_xlim(0, 1)
    ax2.xaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax2.grid(False)
    ax2.set_title("...but most stars that explode are core collapse")
    fig.tight_layout()
    caption(fig, SOURCE_BTS + " Left: classified events weighted by 1/(classification efficiency "
            "at their magnitude). Right: Li et al.'s grouping (II includes IIb and IIn; Ib/c "
            "includes Ic-BL, Ibn, Icn); BTS volumetric shares count supernovae brighter than "
            "M = −14, weighted by 1/V_max. LOSS (Lick Observatory Supernova Search) values from "
            "Li et al. 2011, Table 7.")
    return save(fig, _out("04_magnitude_vs_volume"))


# ---------------------------------------------------------------------------
# 5. Luminosity functions
# ---------------------------------------------------------------------------

def fig_luminosity_functions(res, df) -> list:
    s = res["rates"]["sample_main"]
    lfs = res["rates"]["lf"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.0))
    bins = np.arange(-23, -11.99, 0.5)
    for fam in ("Ia", "II", "Ibc", "IIn", "SLSN"):
        m = s.loc[s["family"] == fam, "M_peak"]
        h, e = np.histogram(m, bins=bins)
        ax1.step(e[:-1], h / h.sum() / 0.5, where="post", color=FAMILY_COLOR[fam], lw=1.5,
                 label=f"{C.FAMILY_LABEL[fam].split(' (')[0]} ({len(m):,})")
    ax1.set_xlim(-12.5, -22.5)
    ax1.set_xlabel("peak absolute magnitude")
    ax1.set_ylabel("fraction of each class per magnitude")
    ax1.set_title("What the survey finds (magnitude-limited)")
    ax1.legend(loc="upper left")
    offs = {"Ia": -0.06, "II": -0.02, "Ibc": 0.02, "IIn": 0.06, "SLSN": 0.0}
    for fam in ("Ia", "II", "Ibc", "IIn", "SLSN"):
        t = lfs[fam]
        mid = 0.5 * (t["M_lo"] + t["M_hi"]).to_numpy() + offs[fam]
        n = t["n"].to_numpy()
        phi, lo, hi = t["phi"].to_numpy(), t["phi_lo"].to_numpy(), t["phi_hi"].to_numpy()
        use = n >= 2
        ax2.plot(mid[use], phi[use], "-", color=FAMILY_COLOR[fam], lw=1.0, alpha=0.6)
        ax2.errorbar(mid[use], phi[use], yerr=[np.clip(phi - lo, 0, None)[use] * 0.999,
                                              np.clip(hi - phi, 0, None)[use]],
                     fmt="none", ecolor=FAMILY_COLOR[fam], lw=1.0, alpha=0.8)
        solid = use & (n >= 5)
        ax2.plot(mid[solid], phi[solid], "o", color=FAMILY_COLOR[fam], ms=4.5, mec=SURFACE, mew=0.6)
        thin = use & (n < 5)
        ax2.plot(mid[thin], phi[thin], "o", mfc=SURFACE, mec=FAMILY_COLOR[fam], ms=4.5, mew=1.0)
    for fam, (M0, _) in (("Ia", lit.LI_2011_MEAN_M["Ia"]), ("II", lit.LI_2011_MEAN_M["II"]),
                         ("Ibc", lit.LI_2011_MEAN_M["Ibc"])):
        ax2.plot([M0], [6e5], "v", color=FAMILY_COLOR[fam], ms=6)
    ax2.text(-17.3, 1.2e6, "▼ mean of the LOSS volume-limited sample (Li et al. 2011)",
             fontsize=7.5, color=INK_SECONDARY, ha="center")
    ax2.set_yscale("log")
    ax2.set_ylim(0.3, 3e6)
    ax2.set_xlim(-12.5, -22.5)
    ax2.set_xlabel("peak absolute magnitude")
    ax2.set_ylabel("rate per magnitude (Gpc⁻³ yr⁻¹ mag⁻¹)")
    ax2.set_title("What actually explodes (per unit volume)")
    for fam, x, y in (("Ia", -19.75, 5e4), ("II", -15.4, 6e4), ("Ibc", -16.8, 1.5e3),
                      ("IIn", -19.6, 3.5e2), ("SLSN", -21.7, 25)):
        ax2.text(x, y, C.FAMILY_LABEL[fam].split(" (")[0], color=INK, fontsize=8, ha="center")
    ax2.plot([], [], "o", color=INK_SECONDARY, ms=4.5, label="5 or more supernovae in the bin")
    ax2.plot([], [], "o", mfc=SURFACE, mec=INK_SECONDARY, ms=4.5, label="2-4 supernovae")
    ax2.legend(loc="lower left", fontsize=7.5)
    fig.tight_layout()
    caption(fig, SOURCE_BTS + " Right: 1/V_max luminosity functions, June 2019 - October 2025, "
            "half-magnitude bins with 68% bootstrap intervals (bins with a single event omitted); "
            "distances from redshifts (CMB frame) or Cosmicflows-4 below 4,000 km/s; host-galaxy "
            "dust not removed. Each faint bin rests on a few nearby events: uncertain by factors "
            "of several.")
    return save(fig, _out("05_luminosity_functions"))


# ---------------------------------------------------------------------------
# 6. Rates
# ---------------------------------------------------------------------------

def fig_rates(res, df) -> list:
    rr = res["rates"]
    s = rr["sample_main"]
    fig = plt.figure(figsize=(13, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.15, 1])
    ax = fig.add_subplot(gs[0])
    grid = np.arange(-22.5, -13.49, 0.25)
    rng = np.random.default_rng(C.RNG_SEED)
    for lab, mask, col in (("Type Ia", s["family"] == "Ia", BLUE), ("core collapse", s["is_cc"], ORANGE)):
        w, M = s.loc[mask, "w"].to_numpy(), s.loc[mask, "M_peak"].to_numpy()
        cum = np.array([w[M < g].sum() for g in grid])
        boots = []
        for _ in range(300):
            k = rng.integers(0, len(w), len(w))
            boots.append([w[k][M[k] < g].sum() for g in grid])
        lo, hi = np.quantile(boots, [0.16, 0.84], axis=0)
        ax.fill_between(grid, lo, hi, color=col, alpha=0.2, lw=0)
        ax.plot(grid, cum, color=col, lw=1.6)
        xl = -19.6 if col == BLUE else -15.2
        ax.text(xl, np.interp(xl, grid, cum) * (2.2 if col == BLUE else 1.9), lab, ha="center",
                fontsize=8, color=INK)
    ax.set_yscale("log")
    ax.set_xlim(-22.5, -13.0)
    ax.set_ylim(1, 3e5)
    ax.set_xlabel("counting supernovae brighter than M")
    ax.set_ylabel("volumetric rate (Gpc⁻³ yr⁻¹)")
    ax.set_title("Rate of supernovae brighter than M")
    ax.invert_xaxis()

    ax2 = fig.add_subplot(gs[1])
    main = rr["main"]["rates"]
    rep = rr["perley_replication"]
    rows = [
        ("Type Ia", "this work (2019-25)", main["Ia"]["rate"], main["Ia"]["lo"], main["Ia"]["hi"], BLUE, True),
        ("Type Ia", "Perley's method, 2018-20", rep["Ia"]["rate"], rep["Ia"]["lo"], rep["Ia"]["hi"], BLUE, False),
        ("Type Ia", "Perley et al. 2020 (BTS)", lit.PERLEY_2020["rate_Ia"][0],
         lit.PERLEY_2020["rate_Ia"][0] - lit.PERLEY_2020["rate_Ia"][1],
         lit.PERLEY_2020["rate_Ia"][0] + lit.PERLEY_2020["rate_Ia"][1], MUTED, False),
        ("Type Ia", "Frohmaier et al. 2019 (PTF)", lit.FROHMAIER_2019_IA[0],
         lit.FROHMAIER_2019_IA[0] - 1.96 * lit.FROHMAIER_2019_IA[1],
         lit.FROHMAIER_2019_IA[0] + 1.96 * lit.FROHMAIER_2019_IA[1], MUTED, False),
        ("Type Ia", "Desai et al. 2026 (ASAS-SN)", lit.DESAI_2026_IA[0],
         lit.DESAI_2026_IA[0] - 1.96 * lit.DESAI_2026_IA[1],
         lit.DESAI_2026_IA[0] + 1.96 * lit.DESAI_2026_IA[1], MUTED, False),
        ("Core collapse", "this work (2019-25)", main["CC"]["rate"], main["CC"]["lo"], main["CC"]["hi"], ORANGE, True),
        ("Core collapse", "our CC/Ia ratio × PTF Ia rate", rr["cc_anchored"]["PTF (Frohmaier+2019)"]["rate"],
         rr["cc_anchored"]["PTF (Frohmaier+2019)"]["rate"] - 1.96 * rr["cc_anchored"]["PTF (Frohmaier+2019)"]["err68"],
         rr["cc_anchored"]["PTF (Frohmaier+2019)"]["rate"] + 1.96 * rr["cc_anchored"]["PTF (Frohmaier+2019)"]["err68"],
         ORANGE, False),
        ("Core collapse", "Perley's method, 2018-20", rep["CC"]["rate"], rep["CC"]["lo"], rep["CC"]["hi"], ORANGE, False),
        ("Core collapse", "Perley et al. 2020 (BTS)", lit.PERLEY_2020["rate_CC"][0],
         lit.PERLEY_2020["rate_CC"][0] - lit.PERLEY_2020["rate_CC"][2],
         lit.PERLEY_2020["rate_CC"][0] + lit.PERLEY_2020["rate_CC"][1], MUTED, False),
        ("Core collapse", "Frohmaier et al. 2021 (PTF)", lit.FROHMAIER_2021["rate_CC"][0],
         lit.FROHMAIER_2021["rate_CC"][0] - 1.96 * lit.FROHMAIER_2021["rate_CC"][2],
         lit.FROHMAIER_2021["rate_CC"][0] + 1.96 * lit.FROHMAIER_2021["rate_CC"][1], MUTED, False),
        ("Core collapse", "Pessi et al. 2025 (ASAS-SN)", lit.PESSI_2025["rate_CC"][0],
         lit.PESSI_2025["rate_CC"][0] - 1.96 * lit.PESSI_2025["rate_CC"][2],
         lit.PESSI_2025["rate_CC"][0] + 1.96 * lit.PESSI_2025["rate_CC"][1], MUTED, False),
    ]
    y = 0
    yt, yl = [], []
    for k, (grp, lab, v, lo, hi, col, ours) in enumerate(rows):
        if k and rows[k - 1][0] != grp:
            y += 0.8
        ax2.plot([lo / 1e4, hi / 1e4], [y, y], color=col, lw=2.2 if ours else 1.4)
        ax2.plot([v / 1e4], [y], "o", color=col, ms=6 if ours else 4.5, mec=SURFACE, mew=0.8)
        yt.append(y)
        yl.append(lab)
        y += 1
    ax2.set_yticks(yt)
    ax2.set_yticklabels(yl, fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xscale("log")
    ax2.set_xlim(1.5, 20)
    ax2.set_xticks([2, 3, 5, 10, 20])
    ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%g"))
    ax2.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax2.set_xlabel("rate (10⁴ Gpc⁻³ yr⁻¹), 95% intervals")
    ax2.set_title("Against other surveys")
    ax2.text(19, -0.55, "Type Ia (M < −16.5)", fontsize=8, color=INK, ha="right", va="center")
    ax2.text(19, yt[5] - 0.55, "core collapse (M < −14)", fontsize=8, color=INK, ha="right",
             va="center")

    ax3 = fig.add_subplot(gs[2])
    sc = rr["seasonal_check"]
    x = np.arange(len(sc))
    for key, col, off in (("Ia", BLUE, -0.12), ("CC", ORANGE, 0.12)):
        v = sc[f"rate_{key}"].to_numpy()
        lo, hi = sc[f"{key}_lo68"].to_numpy(), sc[f"{key}_hi68"].to_numpy()
        ref = np.average(v[1:], weights=1 / ((hi[1:] - lo[1:]) / 2) ** 2)
        ax3.errorbar(x + off, v / ref, yerr=[(v - lo) / ref, (hi - v) / ref], fmt="o", color=col,
                     ms=4.5, capsize=0, lw=1.2, label="Type Ia" if key == "Ia" else "core collapse")
    ax3.axhline(1.0, color=AXIS, lw=0.8)
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"{s_[:4]}-{int(s_[2:4]) + 1:02d}" for s_ in sc["season"]], rotation=45, fontsize=7.5)
    ax3.set_ylabel("rate / mean of June 2019 onward")
    ax3.set_title("Is the survey stable? Season by season")
    ax3.set_ylim(0, 2.0)
    ax3.legend(loc="upper right")
    annotate(ax3, "first season:\nfewer found", (0 - 0.12, sc["rate_Ia"].iat[0] /
             np.average(sc["rate_Ia"][1:], weights=1 / ((sc["Ia_hi68"][1:] - sc["Ia_lo68"][1:]) / 2) ** 2)),
             (0.6, 0.35))
    fig.tight_layout()
    caption(fig, SOURCE_BTS + " Rates use Perley et al.'s active sky area (35% of the sky), our "
            "own season-by-season recovery and classification efficiencies, dust losses from the "
            "SFD map, and comoving volumes with cosmological time dilation (Perley et al. used "
            "Euclidean volumes, which lowers their Type Ia rate by ~13%). The last season is "
            "five months long.")
    return save(fig, _out("06_rates"))


# ---------------------------------------------------------------------------
# 7. Host galaxies on the SFR-mass plane
# ---------------------------------------------------------------------------

def fig_hosts(res, df) -> list:
    ho = res["hosts"]
    sel = ho["plane_sample"]
    gal = ho["background_sample"]
    fig = plt.figure(figsize=(13, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15])
    xe, ye = np.linspace(8.5, 12, 71), np.linspace(-3.5, 1.8, 81)
    H2, _, _ = np.histogram2d(gal["lgm_tot_p50"], gal["sfr_tot_p50"], bins=[xe, ye])
    for k, (lab, mask, col) in enumerate((("Type Ia hosts", sel["family"] == "Ia", BLUE),
                                          ("core-collapse hosts", sel["is_cc"], ORANGE))):
        ax = fig.add_subplot(gs[k])
        ax.pcolormesh(xe, ye, np.log10(H2.T + 1), cmap=mpl.colors.LinearSegmentedColormap.from_list(
            "g", [SURFACE, "#d6d4cc", "#9b998f"]), rasterized=True)
        ax.scatter(sel.loc[mask, "lgm_tot_p50"], sel.loc[mask, "sfr_tot_p50"], s=9, color=col,
                   lw=0.3, edgecolor=SURFACE, rasterized=True)
        mm = np.linspace(8.5, 12, 10)
        ax.plot(mm, mm + C.LOG_SSFR_QUENCHED, color=INK_SECONDARY, lw=0.9)
        ax.text(11.85, 11.85 + C.LOG_SSFR_QUENCHED - 0.25, "sSFR = 10⁻¹¹ / yr", ha="right",
                fontsize=7.5, color=INK_SECONDARY, rotation=33)
        q = ho["quenched"]["-11.0"]["Ia" if k == 0 else "II"]
        ax.set_title(f"{lab} ({int(mask.sum()):,})")
        ax.text(0.04, 0.95, "star-forming\ngalaxies", transform=ax.transAxes, fontsize=8,
                color=INK_SECONDARY, va="top")
        ax.text(0.96, 0.05, "quenched\ngalaxies", transform=ax.transAxes, fontsize=8,
                color=INK_SECONDARY, ha="right")
        ax.set_xlim(8.5, 12)
        ax.set_ylim(-3.5, 1.8)
        ax.set_xlabel("stellar mass, log M★ (M☉)")
        if k == 0:
            ax.set_ylabel("star-formation rate, log SFR (M☉/yr)")
    ax3 = fig.add_subplot(gs[2])
    for key, col, lab, mlab in (("per_mass_Ia", BLUE, "Type Ia", "A·M★ + B·SFR (best fit)"),
                                ("per_mass_CC", ORANGE, "core collapse", "∝ star-formation rate")):
        t = ho[key]
        t = t[(t["n"] > 0) & t["ssfr_mid"].notna()]
        x = t["ssfr_mid"].to_numpy()
        ax3.errorbar(x, t["rel"], yerr=[t["rel"] - t["rel_lo"], t["rel_hi"] - t["rel"]], fmt="o",
                     color=col, ms=5.5, mec=SURFACE, mew=0.7, capsize=0, lw=1.2, label=lab, zorder=3)
        tm = ho[key][ho[key]["ssfr_mid"].notna()]
        ax3.plot(tm["ssfr_mid"], tm["model"], "-", color=col, lw=1.0, alpha=0.7,
                 label=f"model: {mlab}")
    ax3.axhline(1.0, color=AXIS, lw=0.8)
    ax3.set_yscale("log")
    ax3.set_ylim(0.01, 60)
    ax3.set_xlabel("host specific SFR, log (yr⁻¹)")
    ax3.set_ylabel("supernovae per unit stellar mass\n(relative to the average galaxy)")
    ax3.set_title("Per unit of stars: Ia everywhere, CC where stars form")
    ax3.legend(loc="upper left", fontsize=7.5)
    fig.tight_layout()
    qi, qc = ho["quenched"]["-11.0"]["Ia"], ho["quenched"]["-11.0"]["II"]
    caption(fig, "Data: host galaxies matched by directional light radius in SDSS DR17 imaging; "
            "MPA-JHU stellar masses and star-formation rates (Kauffmann et al. 2003; Brinchmann "
            f"et al. 2004). Grey: {len(ho['background_sample']):,} of the SDSS Main Galaxy Sample "
            f"galaxies at 0.015 < z < 0.1. {100 * qi['obs']:.0f}% of Type Ia but {100 * qc['obs']:.0f}% "
            "of Type II hosts lie below the line. Right: hosts per bin divided by the bin's share of the "
            "stellar mass of SDSS galaxies (re-weighted to the hosts' redshifts); 68% binomial "
            "intervals. MPA-JHU SFRs below ~10⁻¹² /yr are poorly constrained.")
    return save(fig, _out("07_host_galaxies"))


def fetch_cutout(ra, dec, scale=0.4, size=256) -> np.ndarray | None:
    """SDSS colour image, cached under data/cutouts."""
    import requests
    from PIL import Image
    d = C.DATA_DIR / "cutouts"
    d.mkdir(exist_ok=True)
    path = d / f"{ra:.5f}_{dec:+.5f}_{scale}_{size}.jpg"
    if not path.exists():
        r = requests.get(C.SKYSERVER_CUTOUT_URL, params={"ra": ra, "dec": dec, "scale": scale,
                                                         "width": size, "height": size}, timeout=60)
        if r.status_code != 200 or len(r.content) < 1000:
            return None
        path.write_bytes(r.content)
    return np.asarray(Image.open(path))


def fig_quenched_hosts(res, df) -> list:
    low = res["hosts"]["lowest_ssfr_cc"]
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.6))
    bpt = {1: "star-forming", 2: "weak SF", 3: "composite", 4: "AGN", 5: "weak LINER", -1: "no lines"}
    for ax, (_, r) in zip(axes.flat, low.iterrows()):
        scale = 0.4
        img = fetch_cutout(r["host_ra"], r["host_dec"], scale=scale)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        if img is None:
            ax.text(0.5, 0.5, "no image", ha="center", transform=ax.transAxes)
            continue
        ax.imshow(img, origin="upper")
        n = img.shape[0]
        dx = -(r["ra"] - r["host_ra"]) * np.cos(np.radians(r["host_dec"])) * 3600 / scale
        dy = -(r["dec"] - r["host_dec"]) * 3600 / scale
        ax.plot(n / 2 + dx, n / 2 + dy, marker="o", ms=14, mfc="none", mec="#ff4f4f", mew=1.3)
        ax.set_title(f"{r['iau_id']} ({r['tns_type'].replace('SN ', '')})", fontsize=9)
        ax.text(0.03, 0.03, f"log sSFR = {r['ssfr']:.1f}\n{bpt.get(int(r['bptclass']), '?')}, "
                f"z = {r['z_sn']:.3f}", transform=ax.transAxes, color="white", fontsize=7.5,
                va="bottom")
    fig.suptitle("Core-collapse supernovae in 'dead' galaxies? The eight lowest-sSFR hosts",
                 x=0.06, ha="left")
    fig.tight_layout(h_pad=2.2)
    caption(fig, "SDSS DR17 colour images, 102″ across, centred on the matched host; red circle: "
            "the supernova. MPA-JHU measures star formation through a 3″ fibre on the galaxy "
            "centre and extrapolates outward from colours; where the centre is an old bulge or a "
            "LINER, a star-forming disc further out can be missed.")
    return save(fig, _out("08_quenched_hosts"))


# ---------------------------------------------------------------------------
# 9. Type II plateaus
# ---------------------------------------------------------------------------

def fig_plateaus(res, df) -> list:
    lc = res["lightcurves"]
    t = lc["table"]
    ok = t[t["plateau_ok"]].copy()
    fig = plt.figure(figsize=(13, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05])
    # examples: spread along the plateau length
    ex = ok.sort_values("fit_t_pt").iloc[np.linspace(2, len(ok) - 3, 3).astype(int)]
    dfi = df.set_index("ztf_id")
    ax0 = fig.add_subplot(gs[0, :2])
    cols = [VIOLET, ORANGE, AQUA]
    for (_, r), col in zip(ex.iterrows(), cols):
        row = dfi.loc[r["ztf_id"]]
        det, _ = LC.load_lc(r["ztf_id"])
        nl = LC.nightly(det)
        z = row["z_cmb"]
        ph = (nl["mjd"] - r["t_exp"]) / (1 + z)
        M = nl["mag"] - row["mu"] - np.where(nl["fid"] == 1, C.A_G_PER_AV, C.A_R_PER_AV) * row["A_V"] \
            + 2.5 * np.log10(1 + z)
        rr_ = nl["fid"] == 2
        ax0.plot(ph[rr_], M[rr_], "o", color=col, ms=3, mec="none")
        tt = np.linspace(20, min(ph.max(), r["fit_t_pt"] + 60), 300)
        ax0.plot(tt, LC.fermi_mag(tt, r["fit_m0"], r["fit_p0"], r["fit_a0"], r["fit_t_pt"], r["fit_w0"]),
                 color=col, lw=1.4)
        ax0.plot([], [], "-o", color=col, ms=3, label=f"{row['iau_id'].replace('SN', 'SN ')}: "
                 f"plateau ends at day {r['fit_t_pt']:.0f}")
    ax0.legend(loc="lower left", fontsize=8)
    ax0.set_xlim(-5, 170)
    ax0.invert_yaxis()
    ax0.set_xlabel("rest-frame days since explosion")
    ax0.set_ylabel("absolute magnitude (ZTF r)")
    ax0.set_title("Plateaus of three Type II supernovae, and the fall onto the cobalt tail")

    ax1 = fig.add_subplot(gs[0, 2])
    ax1.hist(ok["fit_t_pt"], bins=np.arange(50, 140, 7.5), color=ORANGE, alpha=0.8, rwidth=0.92)
    ax1.set_xlabel("plateau end t_PT (days)")
    ax1.set_ylabel("supernovae")
    st = lc["stats"]
    ax1.set_title(f"{st['n_full']} plateaus measured end to end")
    a_mu, a_sd = lit.ANDERSON_2014_OPTD
    top = ax1.get_ylim()[1] * 1.25
    ax1.set_ylim(0, top)
    ax1.errorbar([a_mu + 10], [top * 0.92], xerr=[[a_sd], [a_sd]], fmt="v", color=INK_SECONDARY,
                 ms=6, capsize=0)
    ax1.text(a_mu + 10 - a_sd - 2, top * 0.92, "Anderson et al. 2014 (shifted +10 d:\ntheir 'end' is the "
             "start of the fall)", ha="right", va="center", fontsize=7, color=INK_SECONDARY)

    ax2 = fig.add_subplot(gs[1, 0])
    tail = ok["M_ni_source"] == "tail"
    ax2.scatter(ok.loc[tail, "L50"], ok.loc[tail, "M_ni"], s=18, color=ORANGE, lw=0.3, edgecolor=SURFACE,
                label="nickel from the tail")
    rel = st["ni_relation"]
    xx = np.geomspace(5e40, 1e43, 50)
    ax2.plot(xx, 10 ** (rel["intercept"] + rel["slope"] * np.log10(xx)), color=INK_SECONDARY, lw=1)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("plateau luminosity L₅₀ (erg/s)")
    ax2.set_ylabel("⁵⁶Ni mass (M☉)")
    ax2.set_title("Brighter plateaus, more nickel")
    ax2.axhline(lit.ANDERSON_2019_NI_MEDIAN, color=MUTED, lw=0.8)
    ax2.text(6e40, lit.ANDERSON_2019_NI_MEDIAN * 1.15, "literature median (Anderson 2019)", fontsize=7,
             color=MUTED)

    ax3 = fig.add_subplot(gs[1, 1:])
    for _, r in ok.iterrows():
        ax3.plot([r["M_ej_Rhi"], r["M_ej_Rlo"]], [r["E51_Rhi"], r["E51_Rlo"]], color=GREY_CLOUD, lw=0.8)
    ax3.scatter(ok["M_ej"], ok["E51"], s=20, color=ORANGE, lw=0.3, edgecolor=SURFACE, zorder=3,
                label="R★ = 500 R☉")
    mz = lit.MARTINEZ_2022
    ax3.add_patch(mpl.patches.Rectangle((mz["M_ej_range"][0], mz["E51_range"][0]),
                                        mz["M_ej_range"][1] - mz["M_ej_range"][0],
                                        mz["E51_range"][1] - mz["E51_range"][0], fill=False,
                                        ec=VIOLET, lw=1.2))
    ax3.text(mz["M_ej_range"][1] * 1.06, mz["E51_range"][0] * 1.05, "range of 53 hydrodynamical\n"
             "fits with velocities\n(Martinez et al. 2022)", ha="left", va="bottom", fontsize=7.5,
             color=VIOLET)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlabel("ejected mass M_ej (M☉)")
    ax3.set_ylabel("explosion energy (10⁵¹ erg)")
    ax3.set_title("Weighing the explosions: energy and ejecta from the plateau")
    ax3.text(0.02, 0.96, "grey lines: the same supernova if the star\nwas 1,000 R☉ (left end) or "
             "300 R☉ (right end)", transform=ax3.transAxes, fontsize=7.5, color=INK_SECONDARY, va="top")
    ax3.set_xlim(0.3, 80)
    ax3.set_ylim(0.03, 30)
    fig.tight_layout()
    caption(fig, "Data: ZTF alert photometry via ALeRCE for 1,868 BTS Type II supernovae; "
            f"{st['n_explosion_epoch']:,} have an explosion time (last non-detection within 10 days), "
            f"{st['n_fall_seen']} show the fall off the plateau, {st['n_full']} have both plus g and r "
            "colours at day 50. L₅₀ from Lyman et al. (2014)'s bolometric correction; nickel assumes "
            "full gamma-ray trapping; energies and masses from Goldberg, Bildsten & Paxton (2019)'s "
            "scalings, which leave the progenitor radius free.")
    return save(fig, _out("09_type_II_plateaus"))


# ---------------------------------------------------------------------------
# 10. The red supergiant problem
# ---------------------------------------------------------------------------

def fig_progenitors(res, df) -> list:
    pr = res["progenitors"]
    db = pr["db18_table"].sort_values("logL").reset_index(drop=True)
    s15 = pr["smartt2015_table"].set_index("sn")
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2)

    ax = fig.add_subplot(gs[0, 0])
    for i, r in db.iterrows():
        col = ORANGE if r["kind"] == "detection" else MUTED
        if r["kind"] == "detection":
            ax.errorbar(r["logL"], i, xerr=r["logL_err"], fmt="o", color=col, ms=5, capsize=0)
        else:
            ax.annotate("", xy=(r["logL"] - 0.18, i), xytext=(r["logL"], i),
                        arrowprops={"arrowstyle": "-|>", "color": col, "lw": 1.1})
            ax.plot(r["logL"], i, "|", color=col, ms=8)
        if r["sn"] in s15.index:
            ax.plot(s15.loc[r["sn"], "logL"], i, "x", color=INK_SECONDARY, ms=4.5, mew=1)
        ax.text(3.98, i, r["sn"].replace("SN", "SN "), fontsize=7, va="center", color=INK_SECONDARY)
    ax.axvspan(C.LOGL_HD_LIMIT - 0.02, C.LOGL_HD_LIMIT + 0.02, color=RED, alpha=0.5, lw=0)
    ax.text(C.LOGL_HD_LIMIT + 0.03, len(db) - 1, "brightest red\nsupergiants\n(log L = 5.5)",
            fontsize=7.5, color=RED, va="top")
    ax.axvspan(5.25, 5.5, color=RED, alpha=0.06, lw=0)
    ax.text(5.375, len(db) * 0.42, "no progenitor\nseen here", ha="center", fontsize=7.5, color=RED)
    ax.set_xlim(3.95, 5.85)
    ax.set_yticks([])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("progenitor luminosity, log L/L☉")
    ax.set_title("Stars seen before they exploded (Type II-P/L)")
    ax.plot([], [], "o", color=ORANGE, label="detected (Davies & Beasor 2018)")
    ax.plot([], [], marker=r"$\leftarrow$", ls="", color=MUTED, ms=9, label="upper limit")
    ax.plot([], [], "x", color=INK_SECONDARY, label="Smartt 2015's value (early-M correction)")
    ax.legend(loc="lower right", fontsize=7.5)

    ax2 = fig.add_subplot(gs[0, 1])
    for key, lab, col in (("ml_s15_stars", "Smartt 2015 masses (STARS)", MUTED),
                          ("ml_s15_kepler", "Smartt 2015 masses (KEPLER)", GREY_CLOUD),
                          ("ml_db18", "Davies & Beasor 2018 masses", ORANGE)):
        f = pr[key]
        ax2.plot(f["hi_grid"], f["delta_chi2_hi"], color=col, lw=1.8 if key == "ml_db18" else 1.3,
                 label=f"{lab}: best {f['M_hi']:.1f} M☉")
    for lev, name in ((1, "1σ"), (4, "2σ"), (9, "3σ")):
        ax2.axhline(lev, color=AXIS, lw=0.8)
        ax2.text(39.5, lev + 0.25, name, ha="right", fontsize=7.5, color=MUTED)
    ax2.axvspan(25, 30, color=RED, alpha=0.06, lw=0)
    ax2.text(27.5, 17.5, "where red\nsupergiants\nshould end", ha="center", fontsize=7.5, color=RED)
    ax2.set_xlim(12, 40)
    ax2.set_ylim(0, 20)
    ax2.set_xlabel("heaviest star that explodes as a Type II-P, M_hi (M☉)")
    ax2.set_ylabel("Δχ² (how strongly disfavoured)")
    ax2.set_title("The upper mass limit, fitted the way Smartt et al. did")
    ax2.legend(loc="upper left", fontsize=7.5)

    ax3 = fig.add_subplot(gs[1, 0])
    nb, nbb = pr["nothing_brighter"], pr["nothing_brighter_beasor"]
    bins = np.arange(4.8, 6.0, 0.025)
    ax3.hist(nb["simulated_max"], bins=bins, color=RED, alpha=0.35, density=True,
             label="if red supergiants up to log L = 5.5 exploded")
    ax3.hist(nbb["simulated_max"], bins=bins, histtype="step", color=INK_SECONDARY, density=True,
             lw=1.2, label="the same, with Beasor et al.'s extra 0.2 dex")
    ax3.axvline(nb["brightest_seen"], color=INK, lw=1.6)
    ax3.text(nb["brightest_seen"] - 0.01, ax3.get_ylim()[1] * 0.45,
             f"brightest seen\n(SN 2009hd, {nb['brightest_seen']:.2f})", ha="right", va="top",
             fontsize=8, color=INK)
    ax3.text(0.97, 0.62, f"a sample of {nb['n']} like ours shows nothing\nbrighter than we saw: "
             f"{100 * nb['p']:.0f}% of the time\n({nb['sigma']:.1f}σ; {100 * nbb['p']:.0f}% with the "
             "extra scatter)", transform=ax3.transAxes, ha="right", fontsize=8, color=INK)
    ax3.set_xlabel("brightest progenitor in a simulated sample, log L/L☉")
    ax3.set_ylabel("probability density")
    ax3.set_title("The red supergiant problem as one probability")
    ax3.legend(loc="upper left", fontsize=7.5)
    ax3.set_xlim(4.8, 5.95)

    ax4 = fig.add_subplot(gs[1, 1])
    cvs = res["rates"]["cc_vs_sfr"]
    mmax = np.array([12, 14, 16, 18, 20, 25, 30, 40, 50, 70, 100])
    from analysis import imf_number_per_mass
    kk = [imf_number_per_mass(8, m) for m in mmax]
    ax4.plot(mmax, np.array(kk) * 1e3, color=INK_SECONDARY, lw=1.6,
             label="if every star from 8 M☉ up to M_max explodes (Kroupa IMF)")
    k_lo, k_hi = cvs["k_measured_err68"]
    kj, km = cvs["k_measured_j3"] * 1e3, cvs["k_measured_md14"] * 1e3
    ax4.axhspan(kj / 1.3, km * 1.3, color=ORANGE, alpha=0.12, lw=0)
    ax4.axhspan(kj, km, color=ORANGE, alpha=0.35, lw=0)
    ax4.text(11, kj / 1.3 + 0.08, "± 30%: uncertainty of the star-formation rate", fontsize=7.5,
             color=INK_SECONDARY)
    ax4.text(95, cvs["k_measured_md14"] * 1e3 + 0.1, "measured: core-collapse rate ÷ star-formation\n"
             "rate (Madau & Dickinson 2014 upper edge; project J3 lower)", ha="right", va="bottom",
             fontsize=7.5, color=INK)
    ax4.set_xscale("log")
    ax4.set_xticks([10, 18, 25, 50, 100])
    ax4.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%g"))
    ax4.set_xlabel("M_max: the heaviest star that makes a visible supernova (M☉)")
    ax4.set_ylabel("supernovae per 1,000 M☉ of stars formed")
    ax4.set_title("Can the supernova count weigh the missing stars? Not yet")
    ax4.legend(loc="lower right", fontsize=7.5)
    ax4.set_ylim(4, 13)
    fig.tight_layout()
    caption(fig, "Data: progenitor tables of Smartt (2015) and Davies & Beasor (2018), parsed from "
            "their arXiv sources. Top right: Smartt et al.'s (2009) likelihood with Salpeter IMF; "
            "the lower mass limit fixed by the detections. Bottom left: every supernova's progenitor "
            "drawn from a power law in luminosity with the slope a Salpeter IMF implies (Davies & "
            "Beasor 2020), smeared by its own measurement error; the brightest of each simulated "
            "sample. Bottom right: rates of core-collapse supernovae brighter than M = −14 (this "
            "work, 2019-25) divided by the cosmic star-formation rate at their median redshift; dark "
            "band: the two SFR estimates, light band: their ±30% uncertainty. Supernovae hidden by "
            "dust (19% locally, Mattila et al. 2012; 39% in the infrared SPIRITS survey, Jencson et "
            "al. 2019) or fainter than −14 would raise the measured value.")
    return save(fig, _out("10_red_supergiant_problem"))


ALL = [fig_zoo, fig_gallery, fig_sky, fig_fractions, fig_luminosity_functions, fig_rates,
       fig_hosts, fig_quenched_hosts, fig_plateaus, fig_progenitors]


def draw_all(res, df) -> list:
    out = []
    for f in ALL:
        out += f(res, df)
        print(f"  {f.__name__}")
    return out
