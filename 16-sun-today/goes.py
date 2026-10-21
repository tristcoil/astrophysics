"""GOES X-rays: the week's flux, its flares, and how hot a flare gets.

The X-Ray Sensor (XRS) on NOAA's GOES weather satellites has measured the whole
Sun's soft X-ray output since 1974, in two bands: 1-8 A ("long") and 0.5-4 A
("short"). The quiet Sun barely shines in soft X-rays; a flare heats plasma to 10+
million kelvin and the flux can jump a hundredfold in minutes. A flare's class is
its peak 1-8 A flux: A, B, C, M, X for 1e-8, 1e-7, 1e-6, 1e-5, 1e-4 W/m^2, each
letter ten times the last, with a number within the decade (a C2.5 peaks at
2.5e-6 W/m^2).

Hotter plasma emits relatively more at the shortest wavelengths, so the ratio of
the two bands is a thermometer (Thomas, Starr & Crannell 1985). Given the
temperature, the 1-8 A flux then says how much plasma is emitting: the emission
measure EM = integral of n_e^2 dV (White, Thomas & Schwartz 2005).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import interpolate

import config as C

# --- Flare classes ----------------------------------------------------------------------


def flare_class(flux: float) -> str:
    """GOES class of a peak 1-8 A flux, truncated to one decimal as NOAA does
    (1.376e-6 W/m^2 is a C1.3, not a C1.4)."""
    if not np.isfinite(flux) or flux <= 0:
        return "?"
    letter, base = C.FLARE_CLASSES[0]
    for lt, b in C.FLARE_CLASSES:
        if flux >= b:
            letter, base = lt, b
    return f"{letter}{np.floor(flux / base * 10 + 1e-6) / 10:.1f}"


# --- Loading the flux -----------------------------------------------------------------------


#: SWPC marks eclipse minutes with zero flux, but the minutes either side are
#: already distorted (the Sun is seen through Earth's upper atmosphere). NCEI's
#: science data flag ~11 minutes more on each side; we pad SWPC's gaps to match.
ECLIPSE_MARGIN_MIN = 12


def _load_swpc(path: Path) -> pd.DataFrame:
    """SWPC's operational 1-minute fluxes.

    Dropped: zero values, which mean no data (the satellite in Earth's shadow), not
    zero flux; ``ECLIPSE_MARGIN_MIN`` either side of them; and minutes SWPC flags as
    electron-contaminated, when energetic electrons around the satellite leak into
    the detector and the correction for them dominates the result.
    """
    df = pd.DataFrame(json.loads(path.read_text()))
    df["time"] = pd.to_datetime(df["time_tag"], utc=True)
    cols = {}
    for energy, name in (("0.1-0.8nm", "long"), ("0.05-0.4nm", "short")):
        s = df[df["energy"] == energy].drop_duplicates("time").set_index("time").sort_index()
        s = s.asfreq("1min")
        flux = s["flux"].astype(float)
        dark = ~(flux > 0) & s["flux"].notna()
        pad = dark.rolling(2 * ECLIPSE_MARGIN_MIN + 1, center=True, min_periods=1).max() > 0
        contaminated = s["electron_contaminaton"].fillna(False).astype(bool)
        cols[name] = flux.where(~pad & ~contaminated & (flux > 0))
        cols[f"{name}_all"] = flux.where(~pad & (flux > 0))
    return pd.DataFrame(cols)


def decode_primary_chan(value) -> int:
    """Which detector pair was primary, as the 0-3 index of the response table.

    GOES-R carries two detectors per band, the second (A2, B2) for large flares.
    NCEI packs the choice into one byte: bits 0-1 say A1 (1), A2 (2) or both (3),
    bits 2-3 say B1 (4), B2 (8) or both (12). The response table orders the pairs
    A1+B1, A2+B1, A1+B2, A2+B2. "Both" (a switch within the minute) counts as the
    flare detector.
    """
    v = int(value)
    a = 2 if (v & 3) >= 2 else 1
    b = 2 if ((v & 12) >> 2) >= 2 else 1
    return {(1, 1): 0, (2, 1): 1, (1, 2): 2, (2, 2): 3}[(a, b)]


def _load_ncei(files: list[Path]) -> pd.DataFrame:
    """NCEI's science-quality 1-minute averages. Times are seconds since
    2000-01-01 12:00 UTC ignoring leap seconds, i.e. plain UTC clock time.

    Flags: 0 good, 1 eclipse, 2 bad data, 4 interpolated across a short gap; the
    first two are dropped, interpolated minutes are kept. Minutes where the
    electron-contamination correction exceeds half the observed signal are dropped
    too: the corrected flux is then more correction than Sun.
    """
    import h5py

    parts = []
    for f in files:
        with h5py.File(f, "r") as h:
            t = h["time"][:]
            idx = pd.Timestamp("2000-01-01 12:00:00", tz="UTC") + pd.to_timedelta(t, unit="s")
            d = {
                "long": h["xrsb_flux"][:].astype(float),
                "short": h["xrsa_flux"][:].astype(float),
                "long_flag": h["xrsb_flag"][:].astype(int),
                "short_flag": h["xrsa_flag"][:].astype(int),
                "electron_flag": h["electron_correction_flag"][:].astype(int),
                "primary_chan": [decode_primary_chan(v) for v in h["xrs_primary_chan"][:]],
            }
        parts.append(pd.DataFrame(d, index=idx.round("1min")))
    df = pd.concat(parts).sort_index()
    df = df[~df.index.duplicated()]
    # the short channel before the electron mask: only its pre-flare background is
    # taken from here (see flare_thermal)
    df["short_all"] = df["short"].where(((df["short_flag"] & 3) == 0) & (df["short"] > 0))
    # electron contamination above half the observed signal: bit 16 for XRS-B,
    # bit 8 for XRS-A (NCEI's electron_correction_flag)
    for ch, bit in (("long", 16), ("short", 8)):
        bad = (((df[f"{ch}_flag"] & 3) != 0) | ~(df[ch] > 0)
               | ((df["electron_flag"] & bit) != 0))
        df.loc[bad, ch] = np.nan
    return df.asfreq("1min")


def load_xrs(paths: C.EpochPaths) -> pd.DataFrame:
    """The week's 1-minute X-ray fluxes, science data where it exists.

    Columns ``long``/``short`` (W/m^2) take NCEI's science-quality values where
    they are archived and SWPC's operational values elsewhere (the last day or
    so); ``source`` says which. ``long_swpc``/``long_ncei`` keep both for the
    cross-check, and ``primary_chan`` records which XRS detectors were in use.
    """
    import fetch_data

    m = fetch_data.load_manifest(paths)
    g = m["goes"]
    swpc = _load_swpc(paths.data / g["xrays"])
    ncei_files = [paths.data / v for v in g.get("ncei", {}).values() if v]
    ncei = _load_ncei(ncei_files) if ncei_files else None

    idx = pd.date_range(pd.Timestamp(paths.week_start), pd.Timestamp(paths.epoch), freq="1min")
    out = pd.DataFrame(index=idx)
    out["long_swpc"] = swpc["long"].reindex(idx)
    out["short_swpc"] = swpc["short"].reindex(idx)
    if ncei is not None:
        out["long_ncei"] = ncei["long"].reindex(idx)
        out["short_ncei"] = ncei["short"].reindex(idx)
        out["primary_chan"] = ncei["primary_chan"].reindex(idx)
        covered = ncei["long"].reindex(idx).notna() | ncei["long_flag"].reindex(idx).notna()
    else:
        out["long_ncei"] = out["short_ncei"] = out["primary_chan"] = np.nan
        covered = pd.Series(False, index=idx)
    out["source"] = np.where(covered, "ncei", "swpc")
    out["long"] = np.where(covered, out["long_ncei"], out["long_swpc"])
    out["short"] = np.where(covered, out["short_ncei"], out["short_swpc"])
    short_all_ncei = ncei["short_all"].reindex(idx) if ncei is not None else np.nan
    out["short_all"] = np.where(covered, short_all_ncei, swpc["short_all"].reindex(idx))
    out.attrs["satellite"] = int(g["swpc_satellite"])
    return out


def cross_check(xrs: pd.DataFrame) -> dict:
    """How closely SWPC's operational fluxes track NCEI's science values."""
    out = {}
    for ch in ("long", "short"):
        both = xrs[[f"{ch}_swpc", f"{ch}_ncei"]].dropna()
        if ch == "short":
            both = both[both[f"{ch}_ncei"] > 1e-8]  # below this, electron correction dominates
        if len(both) < 60:
            continue
        ratio = both[f"{ch}_swpc"] / both[f"{ch}_ncei"]
        out[ch] = {
            "n_minutes": int(len(both)),
            "median_ratio": float(ratio.median()),
            "p05_ratio": float(ratio.quantile(0.05)),
            "p95_ratio": float(ratio.quantile(0.95)),
        }
    return out


# --- Flares ----------------------------------------------------------------------------------


def swpc_events(paths: C.EpochPaths) -> pd.DataFrame:
    """SWPC's flare list for the week (events it could not finish, e.g. across an
    eclipse, have no peak and are dropped)."""
    import fetch_data

    m = fetch_data.load_manifest(paths)
    rows = json.loads((paths.data / m["goes"]["flares"]).read_text())
    ev = pd.DataFrame(rows)
    ev = ev.dropna(subset=["begin_time", "max_time", "end_time", "max_class"])
    out = pd.DataFrame({
        "begin": pd.to_datetime(ev["begin_time"], utc=True),
        "peak": pd.to_datetime(ev["max_time"], utc=True),
        "end": pd.to_datetime(ev["end_time"], utc=True),
        "swpc_class": ev["max_class"].values,
        "swpc_peak_flux": ev["max_xrlong"].astype(float).values,
    })
    lo, hi = pd.Timestamp(paths.week_start), pd.Timestamp(paths.epoch)
    return out[(out.peak >= lo) & (out.peak <= hi)].reset_index(drop=True)


def ncei_events(paths: C.EpochPaths) -> pd.DataFrame:
    """NCEI's science flare list (the "flsum" product), for days it has archived.

    Each flare appears as EVENT_START / EVENT_PEAK / EVENT_END rows sharing an id.
    """
    import h5py

    import fetch_data

    m = fetch_data.load_manifest(paths)
    rows = []
    for rel in m["goes"].get("ncei_flares", {}).values():
        if not rel:
            continue
        with h5py.File(paths.data / rel, "r") as h:
            t = pd.Timestamp("2000-01-01 12:00:00", tz="UTC") + pd.to_timedelta(h["time"][:], unit="s")
            status = [s.decode() for s in h["status"][:]]
            cls = [s.decode() for s in h["flare_class"][:]]
            ids = h["flare_id"][:]
        rows += [dict(id=int(i), status=s, time=tt, cls=c) for i, s, tt, c in zip(ids, status, t, cls)]
    if not rows:
        return pd.DataFrame(columns=["begin", "peak", "end", "ncei_class"])
    df = pd.DataFrame(rows)
    out = []
    for fid, g in df.groupby("id"):
        by = {s: g[g.status == s] for s in ("EVENT_START", "EVENT_PEAK", "EVENT_END")}
        if any(v.empty for v in by.values()):
            continue
        out.append({"begin": by["EVENT_START"].time.iloc[0], "peak": by["EVENT_PEAK"].time.iloc[0],
                    "end": by["EVENT_END"].time.iloc[0], "ncei_class": by["EVENT_PEAK"].cls.iloc[0]})
    ev = pd.DataFrame(out).sort_values("peak")
    lo, hi = pd.Timestamp(paths.week_start), pd.Timestamp(paths.epoch)
    return ev[(ev.peak >= lo) & (ev.peak <= hi)].reset_index(drop=True)


def measure_event(xrs: pd.DataFrame, begin, end) -> dict:
    """Background, peak and energy of one flare in the 1-8 A band.

    The background is the median flux in the ``FLARE_BACKGROUND_MIN`` minutes
    before the start; the fluence is the background-subtracted flux integrated from
    start to end (J/m^2), proportional to the energy the flare radiated in 1-8 A.
    """
    pre = xrs["long"].loc[begin - pd.Timedelta(minutes=C.FLARE_BACKGROUND_MIN):
                          begin - pd.Timedelta(minutes=1)]
    bg = float(np.nanmedian(pre)) if pre.notna().any() else np.nan
    seg = xrs["long"].loc[begin:end]
    excess = (seg - bg).clip(lower=0)
    seconds = (seg.index - seg.index[0]).total_seconds().to_numpy()
    good = excess.notna().to_numpy()
    fluence = float(np.trapezoid(excess.to_numpy()[good], seconds[good])) if good.sum() > 1 else np.nan
    return {
        "background": bg,
        "peak_flux": float(np.nanmax(seg)) if seg.notna().any() else np.nan,
        "fluence_J_m2": fluence,
        "missing_minutes": int((~good).sum()),
    }


def flare_list(paths: C.EpochPaths, xrs: pd.DataFrame) -> pd.DataFrame:
    """The week's flares, each re-measured on the flux we use.

    Timing comes from SWPC's real-time event list, which covers the whole week up to
    the epoch. If SWPC's rolling window no longer reaches the epoch's week, NCEI's
    science list takes over. Where NCEI has classified the same flare its class is
    attached, since it comes from the science-quality data.
    """
    ev = swpc_events(paths)
    ncei = ncei_events(paths)
    if ev.empty or ev.begin.min() > pd.Timestamp(paths.week_start) + pd.Timedelta(days=1):
        ev = ncei.assign(swpc_class=None, swpc_peak_flux=np.nan)
        source = "NCEI science flare list"
    else:
        source = "SWPC real-time event list"
        # attach NCEI's class to the same flare (peaks within 3 minutes)
        ncls = []
        for pk in ev.peak:
            near = ncei[(ncei.peak - pk).abs() <= pd.Timedelta(minutes=3)] if len(ncei) else ncei
            ncls.append(near.ncei_class.iloc[0] if len(near) else None)
        ev["ncei_class"] = ncls
    rows = [measure_event(xrs, b, e) for b, e in zip(ev.begin, ev.end)]
    ev = pd.concat([ev.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    ev["class"] = [flare_class(f) for f in ev["peak_flux"]]
    ev.attrs["source"] = source
    return ev


def select_flare(events: pd.DataFrame, xrs: pd.DataFrame, paths: C.EpochPaths) -> dict | None:
    """The week's most energetic flare of class C or above, fully observed."""
    ok = events[(events.peak_flux >= C.FLARE_MIN_CLASS_FLUX) & (events.missing_minutes <= 2)]
    if ok.empty:
        return None
    best = ok.loc[ok.fluence_J_m2.idxmax()]
    return {
        "class": best["class"],
        "swpc_class": best["swpc_class"],
        "ncei_class": best.get("ncei_class"),
        "begin": best.begin.isoformat(),
        "peak": best.peak.isoformat(),
        "end": best.end.isoformat(),
        "peak_flux": float(best.peak_flux),
        "background": float(best.background),
        "fluence_J_m2": float(best.fluence_J_m2),
        "rule": "largest background-subtracted 1-8 A fluence among flares >= C1.0",
        "n_candidates": int(len(ok)),
    }


# --- Temperature and emission measure -----------------------------------------------------


def response_index(satellite: int, primary_chan: int = 0) -> int:
    """Row of the CHIANTI response table for a satellite and detector pair.

    GOES-R satellites (16 on) carry two detectors per band, the second for big
    flares; ``primary_chan`` 0-3 encodes which pair was in use (A1+B1 ... A2+B2).
    """
    if satellite <= 15:
        return satellite - 1
    return 15 + 4 * (satellite - 16) + int(primary_chan)


def temperature_em(short, long, table_path: Path, satellite: int, primary_chan: int = 0,
                   abundance: str = "coronal"):
    """Isothermal temperature (MK) and emission measure (cm^-3) from the two bands.

    The same method and table as SolarSoft's goes_chianti_tem and sunkit-instruments:
    the model flux ratio short/long rises monotonically with temperature, so the
    measured ratio gives T; the 1-8 A flux per unit emission measure at that T then
    gives EM. Fluxes must already be background-subtracted. Points where either band
    is too faint to trust (short < 1e-10 or long < 3e-8 W/m^2) come back as NaN.
    """
    short = np.asarray(short, dtype=float)
    long = np.asarray(long, dtype=float)
    tab = fits.getdata(table_path, 1)
    i = response_index(satellite, primary_chan)
    if abundance == "coronal":
        ratio_model = tab["FSHORT_COR"][i] / tab["FLONG_COR"][i]
        flong_model = tab["FLONG_COR"][i]
    else:
        ratio_model = tab["FSHORT_PHO"][i] / tab["FLONG_PHO"][i]
        flong_model = tab["FLONG_PHO"][i]
    temp_model = tab["TEMP_MK"][i]
    em_scale = 10.0 ** (49.0 - tab["ALOG10EM"][i])  # the table is per 1e49 cm^-3

    t_of_ratio = interpolate.splrep(ratio_model, temp_model, s=0)
    flux_per_em = interpolate.splrep(temp_model, flong_model * em_scale, s=0)

    good = (short >= 1e-10) & (long >= 3e-8) & np.isfinite(short) & np.isfinite(long)
    temp = np.full(short.shape, np.nan)
    em = np.full(short.shape, np.nan)
    ratio = np.where(good, short / np.where(good, long, 1.0), np.nan)
    in_range = good & (ratio >= ratio_model.min()) & (ratio <= ratio_model.max())
    if in_range.any():
        temp[in_range] = interpolate.splev(ratio[in_range], t_of_ratio)
        em[in_range] = long[in_range] / interpolate.splev(temp[in_range], flux_per_em) * 1e49
    return temp, em


def flare_thermal(xrs: pd.DataFrame, sel: dict, table_path: Path) -> pd.DataFrame:
    """Temperature and emission measure through the selected flare, minute by minute.

    Both channels have their pre-flare level subtracted. The short channel's
    pre-flare level is so faint (~5e-9 W/m^2) that electrons dominate it, so it is
    taken from the corrected flux even where flagged; it is 2% of the flare's short
    flux, and ``short_background_sensitivity_MK`` records how far T at the peak moves
    if it is set to zero instead. T and EM are only computed where the short channel
    is clean and the flare dominates the long channel.
    """
    begin = pd.Timestamp(sel["begin"])
    end = pd.Timestamp(sel["end"])
    t0 = begin - pd.Timedelta(minutes=C.FLARE_BACKGROUND_MIN)
    t1 = end + pd.Timedelta(minutes=C.FLARE_AFTER_MIN)
    seg = xrs.loc[t0:t1, ["long", "short", "short_all", "primary_chan"]].copy()
    pre = xrs.loc[t0:begin - pd.Timedelta(minutes=1)]
    bg_long = float(np.nanmedian(pre["long"]))
    bg_short = float(np.nanmedian(pre["short_all"]))
    seg["long_excess"] = seg["long"] - bg_long
    seg["short_excess"] = seg["short"] - bg_short
    chan = seg["primary_chan"].fillna(0).astype(int)

    def solve(short_excess):
        temp = np.full(len(seg), np.nan)
        em = np.full(len(seg), np.nan)
        for pc in np.unique(chan):
            rows = (chan == pc).to_numpy()
            tt, ee = temperature_em(short_excess[rows], seg["long_excess"].to_numpy()[rows],
                                    table_path, xrs.attrs["satellite"], pc)
            temp[rows], em[rows] = tt, ee
        return temp, em

    temp, em = solve(seg["short_excess"].to_numpy())
    temp0, _ = solve(seg["short"].to_numpy())  # no short-channel background at all
    seg["temperature_MK"] = temp
    seg["em_cm3"] = em
    # only trust T and EM while the flare clearly dominates the background
    weak = seg["long_excess"] < 0.1 * np.nanmax(seg["long_excess"])
    seg.loc[weak, ["temperature_MK", "em_cm3"]] = np.nan
    k = int(np.nanargmax(np.where(weak, np.nan, temp)))
    seg.attrs.update(background_long=bg_long, background_short=bg_short,
                     short_background_sensitivity_MK=float(temp0[k] - temp[k]))
    return seg
