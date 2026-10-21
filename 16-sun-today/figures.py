"""Every figure, drawn from the analysis products.

Image figures are PNG only: an SVG of a solar image would just wrap the same
raster in XML, several MB larger. Charts are PNG + SVG, as in the other projects.
"""

from __future__ import annotations

import warnings
from datetime import datetime

import astropy.units as u
import matplotlib.colors as mpl_colors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.visualization import AsinhStretch, ImageNormalize, LogStretch
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from matplotlib.patches import Circle

import config as C
import limb
import plotstyle as ps
import solar_physics as sp

warnings.filterwarnings("ignore")

import sunpy.visualization.colormaps  # noqa: E402,F401  (registers the SDO colour tables)

ps.use_style()
IMAGE_DPI = 150


# --- Displaying solar images --------------------------------------------------------------------


def channel_norm(m, wavelength: int, *, pct_hi: float = 99.9, crop=None) -> ImageNormalize:
    """Display scaling for one channel: log for the EUV lines (whose brightness spans
    four decades from coronal holes to flaring loops), asinh for the UV and visible
    continua.

    Full disc: the ceiling is a high percentile of the disc and, for EUV, the floor
    sits in the faint corona above the limb, so the off-limb glow stays visible
    without lifting the black sky. For the visible continuum the floor sits above the
    telescope's scattered-light halo (a quarter of the disc brightness in AIA 4500 A).
    With ``crop`` (a zoomed map) both limits come from the zoom itself.
    """
    ch = C.CHANNEL_BY_WAVELENGTH.get(wavelength)
    log = ch is not None and ch.stretch == "log"
    if crop is not None:
        d = np.asarray(crop.data, dtype=float)
        vmin, vmax = np.nanpercentile(d, [1.0, 99.7])
        if log:
            vmin = max(vmin, vmax * 1e-3)
            return ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch(100), clip=True)
        return ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.5), clip=True)
    d = np.asarray(m.data, dtype=float)
    r = limb.fast_radius(m)
    on = r < 0.98
    vmax = float(np.nanpercentile(d[on], pct_hi))
    if log:
        ring = (r > 1.15) & (r < 1.25)
        vmin = max(float(np.nanpercentile(d[ring], 20)), vmax * 1e-3)
        return ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch(1000), clip=True)
    if wavelength == 4500:
        centre = float(np.nanmedian(d[r < 0.3]))
        return ImageNormalize(vmin=0.3 * centre, vmax=1.2 * vmax, stretch=AsinhStretch(0.6),
                              clip=True)
    vmin = float(np.nanpercentile(d[on], 0.5)) * 0.4
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.4), clip=True)


def show_map(ax, m, wavelength, *, norm=None, cmap=None):
    """Draw a map in pixel space with the extent in arcseconds (x west, y north)."""
    d = np.asarray(m.data, dtype=float)
    bl, tr = m.bottom_left_coord, m.top_right_coord
    extent = [bl.Tx.to_value(u.arcsec), tr.Tx.to_value(u.arcsec),
              bl.Ty.to_value(u.arcsec), tr.Ty.to_value(u.arcsec)]
    if wavelength == "hmi":
        cmap = cmap or "gray"
        norm = norm or ImageNormalize(vmin=-250, vmax=250, clip=True)
        # colour it by hand: astropy's normalisation does not pass NaNs through as
        # "bad" pixels, so the sky off the disc is painted explicitly
        rgba = plt.get_cmap(cmap)(np.clip(norm(np.nan_to_num(d)), 0, 1))
        rgba[~np.isfinite(d)] = mpl_colors.to_rgba(ps.DARK_PAGE)
        ax.imshow(rgba, origin="lower", extent=extent, interpolation="lanczos")
    else:
        cmap = cmap or C.CHANNEL_BY_WAVELENGTH[wavelength].cmap
        norm = norm or channel_norm(m, wavelength)
        ax.imshow(d, origin="lower", cmap=cmap, norm=norm, extent=extent, interpolation="lanczos")
    ax.set_facecolor(ps.DARK_PAGE)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    return extent


def rgb_composite(maps: dict, channels=(211, 193, 171)) -> np.ndarray:
    """Three channels into red, green and blue, each with its own log scaling."""
    planes = []
    for wl in channels:
        n = channel_norm(maps[wl], wl)
        planes.append(np.clip(n(np.asarray(maps[wl].data, dtype=float)), 0, 1))
    return np.dstack(planes)


def _t_label(tk: float) -> str:
    """A temperature for a panel label, to the precision it deserves: the ion
    formation temperatures are peaks of broad curves, good to ~0.1 dex."""
    if tk >= 3e5:
        mk = tk / 1e6
        return f"{mk:.1f} MK" if mk < 3 else f"{mk:.0f} MK"
    return f"{float(f'{tk:.1g}'):,.0f} K"


def channel_label(wl: int) -> tuple[str, str]:
    ch = C.CHANNEL_BY_WAVELENGTH[wl]
    temps = " & ".join(_t_label(10 ** lt) for lt in ch.log_t)
    return f"{wl} Å", f"{ch.ions} · {temps}\n{ch.region}"


def _panel_text(ax, title: str, sub: str, *, backing: bool = False) -> None:
    """Channel name (bold) and what it shows, in the panel's top-left corner; with
    ``backing``, on a translucent dark card so it reads over bright images."""
    t1 = TextArea(title, textprops={"color": ps.DARK_INK, "fontsize": 11, "fontweight": "bold"})
    t2 = TextArea(sub, textprops={"color": ps.DARK_INK_SECONDARY, "fontsize": 7.4,
                                  "linespacing": 1.35})
    box = VPacker(children=[t1, t2], align="left", pad=0, sep=3)
    anchored = AnchoredOffsetbox(loc="upper left", child=box, pad=0.35, borderpad=0.45,
                                 frameon=backing)
    if backing:
        anchored.patch.set_facecolor(ps.DARK_PAGE)
        anchored.patch.set_alpha(0.62)
        anchored.patch.set_edgecolor("none")
        anchored.patch.set_boxstyle("round,pad=0.25")
    ax.add_artist(anchored)


def image_grid(nrows: int, ncols: int, *, panel_in: float = 3.7, gap_in: float = 0.05,
               header_in: float = 0.95, footer_in: float = 0.45, margin_in: float = 0.1):
    """A dark figure holding a grid of square image panels, sized to fit exactly."""
    w = 2 * margin_in + ncols * panel_in + (ncols - 1) * gap_in
    h = header_in + footer_in + nrows * panel_in + (nrows - 1) * gap_in
    fig = ps.dark_figure(figsize=(w, h))
    axes = []
    for i in range(nrows):
        for j in range(ncols):
            x = margin_in + j * (panel_in + gap_in)
            y = h - header_in - (i + 1) * panel_in - i * gap_in
            axes.append(fig.add_axes([x / w, y / h, panel_in / w, panel_in / h]))
    return fig, axes, (w, h)


def _earth_marker(ax, m, x, y) -> None:
    """Earth to scale: 12,742 km across, at the Sun's distance."""
    d_arcsec = (12742e3 / m.dsun.to_value(u.m)) * 206264.806
    ax.add_patch(Circle((x, y), d_arcsec / 2, color=ps.DARK_INK, lw=0))
    ax.text(x, y - 40, "Earth\nto scale", ha="center", va="top", fontsize=6.5,
            color=ps.DARK_INK_SECONDARY, linespacing=1.1)


# --- Figure 1: the Sun today --------------------------------------------------------------------


def fig_sun_today(products, results, out) -> list:
    maps, hmi = products["maps"], products["hmi"]
    epoch = datetime.fromisoformat(results["epoch"])
    order = ["hmi", 4500, 1700, 1600, 304, 131, 171, 193, 211, 335, 94, "rgb"]
    fig, axes, (w, h) = image_grid(3, 4, header_in=1.05, footer_in=0.5)
    lim = 1150
    for ax, key in zip(axes, order):
        if key == "hmi":
            show_map(ax, hmi, "hmi")
            _panel_text(ax, "Magnetic field", "HMI, line of sight · white: field\n"
                        "pointing toward us, black: away")
            for tr in _regions_at(products, results, epoch):
                ax.text(tr["x"], tr["y"] + 70, str(tr["noaa"]), color=ps.DARK_INK, fontsize=7,
                        ha="center", va="bottom")
        elif key == "rgb":
            ext = show_map(ax, maps[171], 171)
            ax.images[-1].remove()
            ax.imshow(rgb_composite(maps), origin="lower", extent=ext, interpolation="lanczos")
            _panel_text(ax, "211 + 193 + 171 Å", "as red, green, blue: 2, 1.6\n"
                        "and 0.6 MK gas in one picture")
        else:
            show_map(ax, maps[key], key)
            title, sub = channel_label(key)
            if key == 4500:
                sub += "\n(streaks: stray light in the telescope)"
            _panel_text(ax, title, sub)
            if key == 171:
                _earth_marker(ax, maps[key], 960, -1010)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    fig.text(0.1 / w, 1 - 0.12 / h, f"The Sun on {epoch:%-d %B %Y}, {epoch:%H:%M} UTC, in "
             "twelve views", color=ps.DARK_INK, fontsize=17, fontweight="bold", va="top")
    fig.text(0.1 / w, 1 - 0.55 / h, "Ordered by the temperature of the gas each channel sees: "
             "from the 5,800 K surface (top) up through the chromosphere to the million-degree "
             "corona and 6-million-degree\nflare plasma (bottom). Solar north is up; west, the "
             "limb the Sun rotates toward, is right.", color=ps.DARK_INK_SECONDARY,
             fontsize=9.5, va="top", linespacing=1.4)
    fig.text(0.5, 0.2 / h, "SDO/AIA near-real-time synoptic images (1024 x 1024, 2.4\" pixels; "
             "JSOC, Stanford) and an HMI magnetogram, all within a minute of 20:00 UTC. EUV "
             "channels on a log scale, UV/visible on asinh. Numbers mark NOAA active regions. "
             "Ions and temperatures: Lemen et al. (2012).", color=ps.DARK_MUTED, fontsize=7.5,
             ha="center", va="center")
    return ps.save(fig, out / "01_sun_today", formats=("png",), dpi=IMAGE_DPI)


def _regions_at(products, results, epoch) -> list:
    """Tracked NOAA regions at the epoch, positioned for the snapshot images."""
    out = []
    named = {s["noaa"] for s in results["regions"]["noaa_regions"]}
    for tr in products["tracks"]:
        if tr.get("noaa") not in named:
            continue
        pts = tr["track"]
        last = pts[-1]
        dt = abs((datetime.fromisoformat(last["time"]) - epoch.replace(tzinfo=None)).total_seconds())
        if dt < 4 * 3600:
            out.append({"noaa": tr["noaa"], "x": last["hpc_x"], "y": last["hpc_y"]})
    # one label per region
    seen, uniq = set(), []
    for o in out:
        if o["noaa"] not in seen:
            seen.add(o["noaa"])
            uniq.append(o)
    return uniq


# --- Figure 2: one active region, layer by layer -----------------------------------------------


def fig_layers(products, results, out) -> list:
    from scipy import ndimage

    maps, hmi = products["maps"], products["hmi"]
    epoch = datetime.fromisoformat(results["epoch"])
    regs = _regions_at(products, results, epoch)
    info = {s["noaa"]: s for s in results["regions"]["noaa_regions"]}
    # the biggest region on the disc at the epoch
    target = max(regs, key=lambda r: info[r["noaa"]]["max_area_msh"])
    x0, y0, half = target["x"], target["y"], 190.0
    order = ["hmi", 4500, 1700, 1600, 304, 131, 171, 193, 211, 335, 94, "171+hmi"]
    fig, axes, (w, h) = image_grid(3, 4, panel_in=3.4, header_in=1.0, footer_in=0.45)

    def crop(m):
        bl = SkyCoord((x0 - half) * u.arcsec, (y0 - half) * u.arcsec, frame=m.coordinate_frame)
        trc = SkyCoord((x0 + half) * u.arcsec, (y0 + half) * u.arcsec, frame=m.coordinate_frame)
        return m.submap(bl, top_right=trc)

    hsub = crop(hmi)
    for ax, key in zip(axes, order):
        if key == "hmi":
            show_map(ax, hsub, "hmi", norm=ImageNormalize(vmin=-600, vmax=600, clip=True))
            _panel_text(ax, "Magnetic field", "the root of it all: the two\n"
                        "polarities of a bipolar region", backing=True)
        elif key == "171+hmi":
            sub = crop(maps[171])
            show_map(ax, sub, 171, norm=channel_norm(maps[171], 171, crop=sub))
            hb = hsub.bottom_left_coord, hsub.top_right_coord
            ext = [hb[0].Tx.to_value(u.arcsec), hb[1].Tx.to_value(u.arcsec),
                   hb[0].Ty.to_value(u.arcsec), hb[1].Ty.to_value(u.arcsec)]
            sm = ndimage.gaussian_filter(np.nan_to_num(hsub.data), 2.0)
            for level, color in ((100, ps.DARK_INK), (-100, ps.DARK_BLUE)):
                ax.contour(sm, levels=[level], colors=[color], linewidths=1.0, extent=ext,
                           origin="lower", linestyles="solid")
            _panel_text(ax, "171 Å + field", "loops arc from one polarity\n(white: toward "
                        "us) to the other (blue: away)", backing=True)
        else:
            sub = crop(maps[key])
            show_map(ax, sub, key, norm=channel_norm(maps[key], key, crop=sub))
            title, sub_text = channel_label(key)
            _panel_text(ax, title, sub_text, backing=True)
        ax.set_xlim(x0 - half, x0 + half)
        ax.set_ylim(y0 - half, y0 + half)
    s = info[target["noaa"]]
    mm = 2 * half * 0.7255 * (results["snapshot"]["dsun_au"])
    fig.text(0.1 / w, 1 - 0.12 / h, f"Active region NOAA {target['noaa']}, from the surface up",
             color=ps.DARK_INK, fontsize=16, fontweight="bold", va="top")
    fig.text(0.1 / w, 1 - 0.52 / h, f"{epoch:%-d %B %Y %H:%M} UTC. The same {2 * half:.0f}\" "
             f"({mm:.0f} Mm, 22 Earths) square in every panel, at {abs(s['lat']):.0f}°"
             f"{'N' if s['lat'] >= 0 else 'S'}. The field is concentrated in a few spots and "
             "plage at the surface;\nabove it, hot gas is held in loops that reach ever higher "
             "into the corona.", color=ps.DARK_INK_SECONDARY, fontsize=9.5, va="top",
             linespacing=1.4)
    fig.text(0.5, 0.2 / h, "SDO/AIA synoptic images (2.4\" pixels) and HMI magnetogram, as in "
             "figure 1; each panel is scaled to its own range (log for EUV, asinh for UV and "
             "visible). Contours: magnetogram smoothed over 4\", at +-100 G.",
             color=ps.DARK_MUTED, fontsize=7.5, ha="center", va="center")
    return ps.save(fig, out / "02_active_region_layers", formats=("png",), dpi=IMAGE_DPI)


# --- Figure 3: a week of X-rays ------------------------------------------------------------------


def _class_labels(ax, lo: float, hi: float) -> None:
    """GOES class letters down the right edge, one per decade of flux."""
    import matplotlib.transforms as mtransforms

    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for letter, base in C.FLARE_CLASSES:
        if lo <= base * 3.16 <= hi:
            ax.text(1.008, base * 3.16, letter, transform=trans, va="center", ha="left",
                    fontsize=9, fontweight="bold", color=ps.INK_SECONDARY)


def fig_xray_week(products, results, out) -> list:
    xrs, ev = products["xrs"], products["events"]
    xr = results["xrays"]
    sel = results["flare"]["selection"]
    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    ax.plot(xrs.index, xrs["long"], color=ps.BLUE, lw=1.1, zorder=3)
    lo, hi = 5e-8, 1e-5
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ax.set_xlim(xrs.index[0], xrs.index[-1])
    _class_labels(ax, lo, hi)
    # real-time stretch (no science-quality data yet)
    rt = xrs.index[xrs["source"] == "swpc"]
    if len(rt):
        ax.axvspan(rt[0], rt[-1], color=ps.C_CONTEXT, alpha=0.35, lw=0, zorder=0)
        ax.text(rt[0] + (rt[-1] - rt[0]) / 2, 7e-6, "real-time data\n(SWPC)", ha="center",
                va="top", fontsize=7.5, color=ps.INK_SECONDARY)
    # flares of class C and above, labelled; alternate heights where they crowd
    big = ev[ev.peak_flux >= 1e-6].sort_values("peak")
    last_t, bump = None, 0
    for _, e in big.iterrows():
        near = last_t is not None and (e.peak - last_t) < pd.Timedelta(hours=5)
        bump = (bump + 1) % 2 if near else 0
        ax.annotate(e["class"], (e.peak, e.peak_flux), xytext=(0, 5 + 10 * bump),
                    textcoords="offset points", ha="center", va="bottom", fontsize=6.8,
                    color=ps.INK_SECONDARY)
        last_t = e.peak
    pk = pd.Timestamp(sel["peak"])
    reg = results["flare"]["location"].get("noaa") or {}
    ps.annotate(ax, f"flare of the week: {sel['class']},\nNOAA {reg.get('region', '?')} "
                "on the west limb", (pk, sel["peak_flux"] * 1.35),
                (pk - pd.Timedelta(hours=40), 6.5e-6), fontsize=8, color=ps.INK, va="center")
    ax.set_ylabel("1–8 Å X-ray flux at Earth (W/m²)")
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%-d %b"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[6, 12, 18]))
    cls = xr["counts_by_class"]
    ax.set_title(f"A week of solar X-rays: {xr['n_events']} flares, {cls['C']} of them C-class, "
                 f"none bigger than {xr['largest']['class']}", loc="left")
    ps.caption(fig, f"GOES-{xr['satellite']} X-Ray Sensor, 1–8 Å channel (the one flare classes "
               "are defined on), 1-minute averages: NOAA NCEI "
               "science-quality data where archived, SWPC real-time data for the last day. "
               "Flare list: SWPC; classes re-measured on the plotted flux. Gaps: the satellite "
               "in Earth's shadow (daily near 09 UTC around the equinox) and minutes when "
               "energetic electrons dominate the detectors. Background level between flares: "
               f"{xr['background_class']}.")
    fig.tight_layout()
    return ps.save(fig, out / "03_xray_week")


# --- Figure 4: the flare, in X-rays and EUV ---------------------------------------------------------


FLARE_LINE_CHANNELS = (131, 94, 335, 211, 193, 171)
#: One colour per entity across the whole flare figure: the two GOES channels take the
#: first two slots, the six AIA channels (hot to cool) the next six, in fixed order.
FLARE_COLORS = {"long": ps.BLUE, "short": ps.ORANGE, 131: ps.AQUA, 94: ps.YELLOW,
                335: ps.MAGENTA, 211: ps.GREEN, 193: ps.VIOLET, 171: ps.RED}


def fig_flare(products, results, out) -> list:
    fp = products["flare"]
    fl = results["flare"]
    th, curves = fp["thermal"], fp["curves"]
    sel, en = fl["selection"], fl["energetics"]
    begin, end, peak = (pd.Timestamp(sel[k]) for k in ("begin", "end", "peak"))
    t_end = pd.Timestamp(fl["episode_end"])
    xrs = products["xrs"]
    t0, t1 = begin - pd.Timedelta(minutes=15), max(df.index[-1] for df in curves.values())

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(4, 2, width_ratios=[2.25, 1], height_ratios=[1.15, 0.75, 0.75, 1.5],
                          hspace=0.1, wspace=0.22, left=0.07, right=0.97, top=0.865, bottom=0.08)
    ax_x = fig.add_subplot(gs[0, 0])
    ax_t = fig.add_subplot(gs[1, 0], sharex=ax_x)
    ax_em = fig.add_subplot(gs[2, 0], sharex=ax_x)
    ax_lc = fig.add_subplot(gs[3, 0], sharex=ax_x)
    ax_seq = fig.add_subplot(gs[:, 1])

    seg = xrs.loc[t0:t1]
    ax_x.plot(seg.index, seg["long"], color=FLARE_COLORS["long"], lw=1.3)
    ax_x.plot(seg.index, seg["short_all"], color=FLARE_COLORS["short"], lw=1.0)
    ax_x.set_yscale("log")
    ax_x.set_ylim(1e-9, 1e-5)
    _class_labels(ax_x, 1e-9, 1e-5)
    ax_x.text(peak + pd.Timedelta(minutes=4), 2.6e-6, "1–8 Å", color=ps.BLUE, fontsize=8)
    ax_x.text(peak + pd.Timedelta(minutes=4), 2.3e-7, "0.5–4 Å", color=ps.ORANGE, fontsize=8)
    ax_x.set_ylabel("flux (W/m²)")

    ax_t.plot(th.index, th["temperature_MK"], color=ps.INK, lw=1.4)
    ax_t.set_ylabel("temperature\n(MK)")
    ax_t.set_ylim(6, 12)
    ax_em.plot(th.index, th["em_cm3"] / 1e48, color=ps.INK, lw=1.4)
    ax_em.set_ylabel("emission measure\n(10⁴⁸ cm⁻³)")
    ax_em.set_ylim(0, 1.3)
    tt, te = pd.Timestamp(en["time_t_max"]), pd.Timestamp(en["time_em_max"])
    ax_t.plot([tt], [en["t_max_MK"]], "o", color=ps.INK, ms=5)
    ps.annotate(ax_t, f"hottest, {en['t_max_MK']:.1f} MK, "
                f"{abs(fl['t_peak_before_flux_peak_min']):.0f} min before the X-ray peak",
                (tt, en["t_max_MK"]), (18, -2), textcoords="offset points", va="center")
    ax_em.plot([te], [en["em_cm3"] / 1e48], "o", color=ps.INK, ms=5)
    ps.annotate(ax_em, f"most plasma, {en['em_cm3'] / 1e48:.2f}×10⁴⁸ cm⁻³, "
                f"{fl['em_peak_after_flux_peak_min']:.0f} min after it",
                (te, en["em_cm3"] / 1e48), (18, -8), textcoords="offset points", va="center")

    for wl in FLARE_LINE_CHANNELS:
        df = curves[wl]
        ep = df["excess"].loc[:t_end]
        ch = C.CHANNEL_BY_WAVELENGTH[wl]
        lt = ch.log_t[1] if wl == 131 else ch.log_t[0]
        ax_lc.plot(df.index, df["excess"] / ep.max(), color=FLARE_COLORS[wl], lw=1.4,
                   label=f"{wl} Å  {_t_label(10 ** lt)}")
    ax_lc.legend(loc="upper right", ncol=2, fontsize=7.5, frameon=True, facecolor=ps.SURFACE,
                 edgecolor="none", framealpha=0.9, title="AIA channel, hot to cool",
                 title_fontsize=7.5)
    ax_lc.set_ylim(-0.25, 1.75)
    ax_lc.set_ylabel("brightening in the flare\nloops (peak = 1)")
    for ax in (ax_x, ax_t, ax_em, ax_lc):
        ax.axvline(peak, color=ps.MUTED, lw=0.8, ls=":")
        ax.axvspan(t_end, t1, color=ps.C_CONTEXT, alpha=0.4, lw=0)
    ax_x.text(t_end + (t1 - t_end) / 2, 5e-6, "renewed heating in the same region\n"
              "(after 17:05; not used)", ha="center", va="top", fontsize=7.5,
              color=ps.INK_SECONDARY)
    ax_lc.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax_lc.set_xlim(t0, t1)
    for ax in (ax_x, ax_t, ax_em):
        plt.setp(ax.get_xticklabels(), visible=False)
    ax_lc.set_xlabel(f"{peak:%-d %B %Y}, UTC")

    # the cooling sequence: when each channel peaks against the temperature it sees
    for r in fl["cooling_sequence"]:
        if r["kind"] == "chromospheric":
            continue
        tmk = 10 ** r["log_t"] / 1e6
        x = max(r["delay_min"], 0.5)
        col = FLARE_COLORS[r["wavelength"]]
        if r["lower_limit"]:
            ax_seq.annotate("", xy=(x * 1.9, tmk), xytext=(x, tmk),
                            arrowprops={"arrowstyle": "-|>", "color": col, "lw": 1.3})
            ax_seq.plot([x], [tmk], "|", color=col, ms=10, mew=1.6)
            label = f"{r['wavelength']} Å  > {r['delay_min']:.0f} min"
        else:
            lo_, hi_ = r["plateau_min"]
            ax_seq.errorbar([x], [tmk], xerr=[[max(x - max(lo_, 0.5), 0)], [max(hi_ - x, 0)]],
                            fmt="o", color=col, ms=7, capsize=0, elinewidth=1.2)
            label = f"{r['wavelength']} Å  {r['delay_min']:.0f} min"
        ax_seq.text(x * 1.12, tmk * 1.07, label, fontsize=7.5, color=ps.INK_SECONDARY,
                    va="bottom")
    tau = en["radiative_cooling_time_min"]
    ax_seq.axvline(tau, color=ps.MUTED, lw=1.1, ls="--")
    ax_seq.text(tau * 0.93, 13.5, f"radiative cooling\ntime at the flare's\ndensity: "
                f"{tau:.0f} min", color=ps.INK_SECONDARY, fontsize=7.5, va="top", ha="right")
    ax_seq.set_xscale("log")
    ax_seq.set_yscale("log")
    ax_seq.set_xlim(0.8, 600)
    ax_seq.set_ylim(0.45, 16)
    ax_seq.set_yticks([0.5, 1, 2, 5, 10])
    ax_seq.set_yticklabels(["0.5", "1", "2", "5", "10"])
    ax_seq.set_xticks([1, 3, 10, 30, 100, 300])
    ax_seq.set_xticklabels(["1", "3", "10", "30", "100", "300"])
    ax_seq.set_xlabel("minutes after the X-ray peak")
    ax_seq.set_ylabel("temperature the channel sees (MK)")
    ax_seq.set_title("Hotter channels peak first", loc="left", fontsize=10.5)

    reg = fl["location"].get("noaa") or {}
    fig.suptitle(f"The flare of the week: {sel['class']} from NOAA {reg.get('region', '?')}, "
                 f"{peak:%-d %B %H:%M} UTC", x=0.07, ha="left", y=0.985, fontsize=13)
    fig.text(0.07, 0.95, "GOES sees the whole Sun in X-rays (top three panels); AIA sees the "
             f"flare loops themselves (bottom). The flare heats plasma to {en['t_max_MK']:.1f} "
             "million K within minutes; then the loops cool for hours,\nlighting up one AIA "
             "channel after another from hot to cool. Thermal energy at the peak: "
             f"{en['thermal_energy_erg'] / 1e29:.1f}×10²⁹ erg, from "
             f"{en['density_cm3'] / 1e10:.1f}×10¹⁰ electrons/cm³ in an "
             f"{en['volume_cm3'] ** (1 / 3) / 1e8:.0f}-Mm cube of hot loops.",
             fontsize=9, color=ps.INK_SECONDARY, va="top", linespacing=1.45)
    ps.caption(fig, "Temperature and emission measure: isothermal fit to the ratio of the two "
               "GOES-18 channels (NCEI science data, pre-flare level subtracted; CHIANTI response "
               "of White et al. 2005, coronal abundances), shown where the short channel is "
               "free of electron contamination. Light curves: AIA synoptic cut-outs, summed "
               "over the hot loops (131/94 Å brightening above 25% of peak), pre-flare image "
               "subtracted. Error bars: the span within 10% of each peak; arrows: still "
               "rising when the episode ended.", y=0.03)
    return ps.save(fig, out / "04_flare")


def fig_flare_images(products, results, out) -> list:
    import fetch_data
    import flare as flare_mod

    paths = products["epoch_paths"]
    manifest = fetch_data.load_manifest(paths)
    fl = results["flare"]
    loc = fl["location"]
    peak = pd.Timestamp(fl["selection"]["peak"])
    rows = (131, 94, 335, 211, 171)
    offsets = (-16, 0, 11, 41, 80, 104)  # minutes from the X-ray peak
    x0, y0, half = loc["hpc_x_arcsec"] - 25, loc["hpc_y_arcsec"] + 5, 95.0
    fig, axes, (w, h) = image_grid(len(rows), len(offsets), panel_in=1.9, gap_in=0.04,
                                   header_in=1.25, footer_in=0.45, margin_in=0.55)
    mask = products["flare"]["mask"]
    for i, wl in enumerate(rows):
        maps = flare_mod.load_sequence(paths, manifest, wl)
        times = pd.DatetimeIndex([pd.Timestamp(m.date.isot, tz="UTC") for m in maps])
        allv = np.concatenate([np.asarray(m.data).ravel() for m in maps])
        lo_, hi_ = np.nanpercentile(allv, [2, 99.9])
        norm = ImageNormalize(vmin=max(lo_, hi_ * 3e-3), vmax=hi_, stretch=LogStretch(100),
                              clip=True)
        for j, off in enumerate(offsets):
            ax = axes[i * len(offsets) + j]
            k = int(np.argmin(np.abs((times - (peak + pd.Timedelta(minutes=off))).total_seconds())))
            m = maps[k]
            show_map(ax, m, wl, norm=norm)
            ax.set_xlim(x0 - half, x0 + half)
            ax.set_ylim(y0 - half, y0 + half)
            if i == 0:
                dt = (times[k] - peak).total_seconds() / 60
                ax.text(0.5, 1.04, f"{times[k]:%H:%M}\n({dt:+.0f} min)", transform=ax.transAxes,
                        ha="center", va="bottom", fontsize=8, color=ps.DARK_INK_SECONDARY)
            if j == 0:
                ch = C.CHANNEL_BY_WAVELENGTH[wl]
                lt = ch.log_t[1] if wl == 131 else ch.log_t[0]
                ax.text(-0.06, 0.5, f"{wl} Å\n{_t_label(10 ** lt)}", transform=ax.transAxes,
                        ha="right", va="center", fontsize=8.5, color=ps.DARK_INK)
                if i == 0:
                    bl = m.bottom_left_coord
                    sc = m.scale.axis1.to_value(u.arcsec / u.pix)
                    ny, nx = mask.shape
                    xs = bl.Tx.to_value(u.arcsec) + (np.arange(nx) + 0.5) * sc
                    ys = bl.Ty.to_value(u.arcsec) + (np.arange(ny) + 0.5) * sc
                    ax.contour(xs, ys, mask.astype(float), levels=[0.5], colors=[ps.DARK_INK],
                               linewidths=0.8)
    fig.text(0.55 / w, 1 - 0.1 / h, "The flare loops cool through AIA's channels",
             color=ps.DARK_INK, fontsize=14, fontweight="bold", va="top")
    fig.text(0.55 / w, 1 - 0.42 / h, "Rows from hot (top) to cool; each row on one scale, so "
             "brightness changes are real. White contour (top left): the hot loops the light "
             "curves are summed over.", color=ps.DARK_INK_SECONDARY, fontsize=8.5, va="top")
    fig.text(0.5, 0.2 / h, f"SDO/AIA synoptic cut-outs, {peak:%-d %B %Y}; {2 * half:.0f}\" "
             "square on the west limb: the bright vertical edge is the limb, the disc lies to "
             "its left, the loops stand above it against the sky.",
             color=ps.DARK_MUTED, fontsize=7.5, ha="center", va="center")
    return ps.save(fig, out / "05_flare_images", formats=("png",), dpi=IMAGE_DPI)


# --- Figure 6: rotation -------------------------------------------------------------------------------


def fig_rotation(products, results, out) -> list:
    rot = results["rotation"]
    bands = rot["bands"]
    lat = np.array([b["lat"] for b in bands])
    om = np.array([b["omega"] for b in bands])
    se = np.array([b["sem"] for b in bands])
    coef = np.array(rot["fits"]["ABC"]["coef"])
    cov = np.array(rot["fits"]["ABC"]["cov"])
    grid = np.linspace(-50, 50, 201)
    s2 = np.sin(np.radians(grid)) ** 2
    x = np.column_stack([np.ones_like(s2), s2, s2**2])
    fit = x @ coef
    band = np.sqrt(np.einsum("ij,jk,ik->i", x, cov, x))
    pub = C.PUBLISHED_ROTATION

    fig, (ax, axr) = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [1.35, 1]})
    ax.fill_between(grid, fit - band, fit + band, color=ps.ORANGE, alpha=0.18, lw=0)
    ax.plot(grid, fit, color=ps.ORANGE, lw=1.8)
    k = pub["komm1993"]
    ax.plot(grid, sp.rotation_law(grid, k["A"], k["B"], k["C"]), color=ps.VIOLET, lw=1.3)
    d = pub["snodgrass1984"]
    ax.plot(grid, sp.rotation_law(grid, d["A"], d["B"], d["C"]), color=ps.AQUA, lw=1.3)
    b = pub["balthasar1986"]
    gs_ = grid[np.abs(grid) <= 35]
    ax.plot(gs_, sp.rotation_law(gs_, b["A"], b["B"], b["C"]), color=ps.YELLOW, lw=1.3)
    ax.errorbar(lat, om, yerr=se, fmt="o", color=ps.INK, ms=4.5, elinewidth=1, capsize=0,
                zorder=5)
    for s in results["regions"]["noaa_regions"]:
        if "omega" in s and s.get("omega_err", 9) < 0.3:
            ax.errorbar([s["lat"]], [s["omega"]], yerr=[s["omega_err"]], fmt="D", ms=5,
                        mfc=ps.SURFACE, mec=ps.BLUE, ecolor=ps.BLUE, elinewidth=1, capsize=0,
                        zorder=6)
            ax.text(s["lat"] + 1.2, s["omega"] + 0.04, str(s["noaa"]), fontsize=6.8,
                    color=ps.BLUE, zorder=6)
    ax.text(48, sp.rotation_law(48, *coef) - 0.12, "this week (fit)", color=ps.ORANGE,
            fontsize=8, ha="right", va="top")
    ps.annotate(ax, "magnetic features,\nKomm et al. 1993",
                (40, sp.rotation_law(40, k["A"], k["B"], k["C"])), (14, 12.85),
                color=ps.VIOLET, fontsize=8, va="center")
    ax.text(-8, 14.1, "surface gas (Doppler),\nSnodgrass 1984", color=ps.AQUA, fontsize=8,
            ha="center", va="top")
    ax.text(-35, sp.rotation_law(35, b["A"], b["B"]) - 0.05, "sunspots,\nBalthasar et al. 1986",
            color=ps.YELLOW, fontsize=8, ha="left", va="top")
    ax.text(28, 14.62, "◇ tracked active regions\n   (NOAA numbers)", color=ps.BLUE,
            fontsize=7.5)
    ax.set_xlim(-52, 52)
    ax.set_ylim(12.6, 14.9)
    ax.set_xlabel("latitude (degrees; south negative)")
    ax.set_ylabel("rotation rate relative to the stars (degrees/day)")
    a, ea = coef[0], rot["fits"]["ABC"]["err"][0]
    ax.set_title(f"This week: the equator turns in {360 / a:.1f} days, latitude 45° in "
                 f"{rot['sidereal_period_45_d']:.1f}", loc="left", fontsize=10.5)

    # right: the difference from the surface-gas (Doppler) rate
    dd = om - sp.rotation_law(lat, d["A"], d["B"], d["C"])
    axr.axhline(0, color=ps.AQUA, lw=1.3)
    axr.plot(grid, sp.rotation_law(grid, k["A"], k["B"], k["C"])
             - sp.rotation_law(grid, d["A"], d["B"], d["C"]), color=ps.VIOLET, lw=1.3)
    axr.plot(grid, fit - sp.rotation_law(grid, d["A"], d["B"], d["C"]), color=ps.ORANGE, lw=1.8)
    axr.errorbar(lat, dd, yerr=se, fmt="o", color=ps.INK, ms=4.5, elinewidth=1, capsize=0)
    cmp_ = rot["comparisons"]["snodgrass1984"]
    axr.text(0, 0.72, f"magnetic network: {cmp_['mean_diff']:+.2f} ± {cmp_['mean_diff_err']:.2f} "
             f"deg/day\n= {cmp_['mean_diff'] / 14.05 * 100:.1f}% faster than the gas around it",
             ha="center", va="top", fontsize=8.5, color=ps.INK)
    axr.set_xlim(-52, 52)
    axr.set_ylim(-0.3, 0.8)
    axr.set_xlabel("latitude (degrees)")
    axr.set_ylabel("rate minus the surface-gas rate (degrees/day)")
    axr.set_title("Magnetic features outrun the surface", loc="left", fontsize=10.5)
    ps.caption(fig, f"{rot['n_frames']} HMI magnetograms (every 3 h for 7 days) remapped to "
               f"Carrington coordinates; {rot['n_pairs']} pairs {rot['pair_separation_h']} h "
               "apart cross-correlated in 5° latitude bands, active regions masked (as Komm et "
               "al. 1993). Error bars: scatter between pairs over the 7 independent days. Fit: "
               f"ω = {coef[0]:.2f} {coef[1]:+.2f} sin²φ {coef[2]:+.2f} sin⁴φ; A = {a:.3f} ± "
               f"{ea:.3f} (stat) ± {rot['A_systematic']:.3f} (pair spacing 6–24 h). Published "
               "laws as compiled by Sudar et al. (2015). Diamonds: individual regions tracked "
               "for 2+ days; each also has its own motion.", y=-0.03)
    fig.tight_layout()
    return ps.save(fig, out / "06_rotation")


# --- Figure 7: the week's field and Hale's law ----------------------------------------------------------


def fig_week_map(products, results, out) -> list:
    syn = products["synoptic"]
    lon, lat, data = syn["lon"], syn["lat"], syn["data"]
    fig, ax = plt.subplots(figsize=(12.5, 5.6))
    cm = plt.get_cmap("gray").copy()
    cm.set_bad(ps.SURFACE)
    ext = [lon[0], lon[-1], lat[0], lat[-1]]
    ax.imshow(np.ma.masked_invalid(data), origin="lower", cmap=cm, vmin=-120, vmax=120,
              extent=ext, aspect="auto", interpolation="nearest")
    ax.grid(False)
    ax.axhline(0, color=ps.MUTED, lw=0.6, ls=":")
    hale = {h["noaa"]: h for h in results["regions"]["hale"]}
    placed = []  # label positions so far, to stagger labels that would collide
    for s in sorted(results["regions"]["noaa_regions"], key=lambda r: r["lon"]):
        if s["noaa"] not in hale or "pos_lon" not in s:
            continue
        h = hale[s["noaa"]]
        lead_pos = h["leading"] > 0
        (lf, bf), (ll, bl) = (((s["neg_lon"], s["neg_lat"]), (s["pos_lon"], s["pos_lat"]))
                              if lead_pos else
                              ((s["pos_lon"], s["pos_lat"]), (s["neg_lon"], s["neg_lat"])))
        lf = lon[0] + (lf - lon[0]) % 360
        ll = lon[0] + (ll - lon[0]) % 360
        if not (lon[0] <= lf <= lon[-1] and lon[0] <= ll <= lon[-1]):
            continue
        color = ps.ORANGE if h["obeys"] else ps.RED
        ax.annotate("", xy=(ll, bl), xytext=(lf, bf),
                    arrowprops={"arrowstyle": "-|>", "color": color, "lw": 2.2,
                                "mutation_scale": 14})
        north = s["lat"] > 0
        lx = (lf + ll) / 2
        ly = (max(bf, bl) + 6) if north else (min(bf, bl) - 6)
        for px, py in placed:
            if abs(px - lx) < 9 and abs(py - ly) < 5:
                ly += 7 if north else -7
        placed.append((lx, ly))
        ax.text(lx, ly, f"{s['noaa']}: {'+' if lead_pos else '−'} leads", ha="center",
                va="bottom" if north else "top", fontsize=7.5, color=ps.INK,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": ps.SURFACE, "alpha": 0.85,
                      "lw": 0})
    ax.set_xlim(lon[0], lon[-1])
    ax.set_ylim(-55, 55)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v % 360:.0f}°"))
    ax.set_xlabel("Carrington longitude (degrees) → west, the direction of rotation")
    ax.set_ylabel("latitude (degrees)")
    hr = results["regions"]
    for y, text, va in ((0.97, "Cycle 25 rule, north: + (white) leads", "top"),
                        (0.03, "Cycle 25 rule, south: − (black) leads", "bottom")):
        ax.text(0.01, y, text, transform=ax.transAxes, fontsize=9, va=va, color=ps.INK,
                fontweight="bold", bbox={"boxstyle": "round,pad=0.3", "facecolor": ps.SURFACE,
                                         "alpha": 0.85, "lw": 0})
    ax.set_title(f"Hale's law on this week's Sun: {hr['hale_obeyed']} of {hr['hale_tested']} "
                 "active regions lead with the polarity Cycle 25 predicts", loc="left")
    ps.caption(fig, "The week's HMI magnetograms merged into one map of the radial field "
               "(each contributes the strip within 25° of its central meridian; white = field "
               "pointing out of the Sun, black = in, saturated at ±120 G). Arrows run from each "
               "region's following polarity to its leading one (flux-weighted centres). In "
               "Solar Cycle 25 the leading polarity is positive in the north and negative in the "
               "south; Cycle 24 was the reverse. Blank: never near central meridian this week.",
               y=-0.04)
    fig.tight_layout()
    return ps.save(fig, out / "07_week_map")


# --- Figures 8 and 9: the limb ------------------------------------------------------------------------------


def fig_limb_darkening(products, results, out) -> list:
    lc = results["limb"]["continuum"]
    mu = np.array(lc["mu"])
    ratio = np.array(lc["ratio"])
    fit = mu >= lc["fit_mu_min"]
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13, 4.4),
                                     gridspec_kw={"width_ratios": [1.2, 0.9, 1.0]})
    mm = np.linspace(0, 1, 200)
    a1.plot(mm, sp.allen_limb_darkening(mm, C.HMI_WAVELENGTH_A), color=ps.VIOLET, lw=1.4)
    a1.plot(mm, sp.eddington_limb_darkening(mm), color=ps.AQUA, lw=1.4)
    a1.plot(mm, 1 + lc["u_fit"] * (mm - 1) + lc["v_fit"] * (mm**2 - 1), color=ps.ORANGE, lw=1.4,
            ls="--")
    a1.plot(mu[fit], ratio[fit], "o", color=ps.INK, ms=3.5)
    a1.plot(mu[~fit], ratio[~fit], "o", mfc="none", mec=ps.MUTED, ms=3.5)
    a1.text(0.02, 0.55, "Allen's Astrophysical\nQuantities, 6173 Å", color=ps.VIOLET,
            fontsize=8, transform=a1.transAxes)
    a1.text(0.3, 0.32, "grey atmosphere\n(Eddington): (2+3μ)/5", color=ps.AQUA, fontsize=8,
            transform=a1.transAxes)
    a1.text(0.55, 0.14, "● HMI, measured", color=ps.INK, fontsize=8, transform=a1.transAxes)
    a1.text(0.55, 0.07, "-- quadratic fit to it", color=ps.ORANGE, fontsize=8,
            transform=a1.transAxes)
    a1.set_xlim(0, 1)
    a1.set_ylim(0.3, 1.02)
    a1.set_xlabel("μ = cos(angle from the vertical)   ← limb · centre →")
    a1.set_ylabel("brightness / brightness at disc centre")
    a1.set_title("The edge of the Sun is darker", loc="left", fontsize=10.5)

    resid = ratio / sp.allen_limb_darkening(mu, C.HMI_WAVELENGTH_A) - 1
    a2.axhline(0, color=ps.VIOLET, lw=1.2)
    a2.plot(mu[fit], 100 * resid[fit], "o", color=ps.INK, ms=3.5)
    a2.plot(mu[~fit], 100 * resid[~fit], "o", mfc="none", mec=ps.MUTED, ms=3.5)
    a2.plot(mm, 100 * (sp.eddington_limb_darkening(mm) / sp.allen_limb_darkening(mm, C.HMI_WAVELENGTH_A) - 1),
            color=ps.AQUA, lw=1.2)
    a2.set_xlim(0, 1)
    a2.set_ylim(-8, 12)
    a2.set_xlabel("μ")
    a2.set_ylabel("difference from Allen's law (%)")
    a2.set_title(f"Measured vs tabulated: {100 * lc['rms_resid_allen']:.1f}% rms", loc="left",
                 fontsize=10.5)

    tau = np.array(lc["tau"])
    a3.plot(tau, lc["t_tau_grey"], color=ps.AQUA, lw=1.4)
    a3.plot(tau, lc["t_tau"], "o", color=ps.INK, ms=3.5)
    a3.set_xscale("log")
    a3.set_xlabel("optical depth τ (how deep we see; τ = μ)")
    a3.set_ylabel("temperature (K)")
    a3.text(0.12, 5700, "grey atmosphere", color=ps.AQUA, fontsize=8)
    a3.text(0.35, 5150, "from this image,\nvia Eddington–Barbier", color=ps.INK, fontsize=8)
    a3.set_title(f"Hotter with depth: {lc['t_at_tau_0p1']:,.0f} K → {lc['t_at_tau_1']:,.0f} K",
                 loc="left", fontsize=10.5)
    ps.caption(fig, f"SDO/HMI continuum image of {lc['date'][:16].replace('T', ' ')} UTC "
               "(4096², definitive; the near-real-time feed has no continuum images). Median "
               "brightness in 60 rings of equal μ; open points (μ < 0.1, the outer 0.5% of the "
               "radius) are blurred by the telescope and excluded from fits. Temperatures assume "
               "local thermodynamic equilibrium and fix T = 5,772 K (the Sun's effective "
               "temperature) at τ = 2/3.", y=-0.03)
    fig.tight_layout()
    return ps.save(fig, out / "08_limb_darkening")


def fig_corona(products, results, out) -> list:
    lb = results["limb"]
    prof, cor = lb["profiles"], lb["corona"]
    lc = lb["continuum"]
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13.5, 4.6),
                                     gridspec_kw={"width_ratios": [1.15, 1.05, 0.9]})
    r = np.array(prof["r"])
    rc = np.sqrt(1 - np.clip(np.array(lc["mu"]), 0, 1) ** 2)
    order = np.argsort(rc)
    a1.plot(rc[order], np.array(lc["ratio"])[order], color=ps.INK, lw=1.6,
            label="visible light (HMI), opaque")
    series = ((1700, ps.BLUE, "1700 Å ultraviolet, opaque"), (304, ps.ORANGE, "304 Å He II"),
              (171, ps.AQUA, "171 Å Fe IX, transparent"), (211, ps.MAGENTA, "211 Å Fe XIV, transparent"))
    for wl, col, lab in series:
        a1.plot(r, np.array(prof[str(wl)]), color=col, lw=1.4, label=lab)
    a1.legend(loc="lower left", fontsize=7.5)
    a1.axvline(1.0, color=ps.MUTED, lw=0.8, ls=":")
    a1.set_yscale("log")
    a1.set_xlim(0, 1.3)
    a1.set_ylim(0.05, 5)
    a1.set_xlabel("distance from disc centre (solar radii)")
    a1.set_ylabel("brightness / disc centre")
    a1.set_title("Dark edge below, bright edge above", loc="left", fontsize=10.5)

    h = np.array(cor["height_mm"])
    cols = {"171": ps.AQUA, "193": ps.YELLOW, "211": ps.MAGENTA, "335": ps.VIOLET}
    for wl, col in cols.items():
        c = cor["channels"][wl]
        prof_q = np.array(c["quiet_profile"])
        base = np.interp(15.0, h, prof_q)
        a2.plot(h, prof_q / base, color=col, lw=1.5)
        hm = np.array(c["quiet_model_heights"])
        a2.plot(hm, np.array(c["quiet_model"]) / base, color=ps.INK, lw=0.9, ls="--")
        lift = {"171": 1.15, "193": 1.15, "211": 0.72, "335": 1.45}[wl]
        a2.text(h[-12], (prof_q / base)[-12] * lift, f"{wl} Å", color=col, fontsize=7.5)
    a2.axvspan(*cor["fit_range_mm"], color=ps.C_CONTEXT, alpha=0.35, lw=0)
    a2.set_yscale("log")
    a2.set_xlim(0, 200)
    a2.set_ylim(0.01, 2)
    a2.set_xlabel("height above the visible surface (Mm)")
    a2.set_ylabel("brightness / brightness at 15 Mm")
    a2.set_title("Hotter gas reaches higher", loc="left", fontsize=10.5)
    a2.text(0.98, 0.95, "quiet corona, median of\n"
            f"{len(cor['quiet_sectors'])} sectors round the limb\n-- hydrostatic fit (shaded range)",
            transform=a2.transAxes, ha="right", va="top", fontsize=7.5, color=ps.INK_SECONDARY)

    for wl, col in cols.items():
        c = cor["channels"][wl]
        x, y = c["t_formation_MK"], c["t_quiet_median_MK"]
        a3.errorbar([x], [y], yerr=[[y - c["t_quiet_p16_MK"]], [c["t_quiet_p84_MK"] - y]],
                    fmt="o", color=col, ms=7, elinewidth=1.2, capsize=0)
        a3.text(x * 1.06, y, f"{wl} Å · {c['ion']}", color=col, fontsize=7.5, va="center")
    tt = np.array([0.4, 3.2])
    a3.plot(tt, tt, color=ps.MUTED, lw=0.9, ls=":")
    a3.text(2.3, 2.45, "equal", color=ps.MUTED, fontsize=7.5, rotation=38)
    a3.set_xscale("log")
    a3.set_yscale("log")
    a3.set_xlim(0.45, 3.2)
    a3.set_ylim(0.45, 3.2)
    for ax in (a3,):
        ax.set_xticks([0.5, 1, 2, 3])
        ax.set_xticklabels(["0.5", "1", "2", "3"])
        ax.set_yticks([0.5, 1, 2, 3])
        ax.set_yticklabels(["0.5", "1", "2", "3"])
    from matplotlib.ticker import NullFormatter
    a3.xaxis.set_minor_formatter(NullFormatter())
    a3.yaxis.set_minor_formatter(NullFormatter())
    a3.set_xlabel("temperature the ion forms at (MK)")
    a3.set_ylabel("temperature from its scale height (MK)")
    a3.set_title("A thermometer from gravity", loc="left", fontsize=10.5)
    ps.caption(fig, "Left: median brightness in rings (AIA snapshot; the visible-light curve is "
               "HMI's continuum). The photosphere is opaque, so its edge shows cooler, higher "
               "layers; the corona is transparent, so its edge adds up the longest path through "
               "it. Middle and right: above the limb each channel's brightness falls as "
               "density² of an isothermal atmosphere in hydrostatic balance in the Sun's 1/r² "
               "gravity (fitted 12–80 Mm, integrated along the line of sight). Error bars: 16–84% "
               "range over quiet sectors; active regions and coronal holes excluded. Formation "
               "temperatures: Lemen et al. (2012).", y=-0.03)
    fig.tight_layout()
    return ps.save(fig, out / "09_corona")


FIGURES = (fig_sun_today, fig_layers, fig_xray_week, fig_flare, fig_flare_images, fig_rotation,
           fig_week_map, fig_limb_darkening, fig_corona)


def make_all(products, results, out) -> list:
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for f in FIGURES:
        written += f(products, results, out)
    return written
