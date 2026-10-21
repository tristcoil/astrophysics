"""Did any runaway start inside a supernova remnant we can still see?

A remnant -- the expanding shell of debris and swept-up gas -- stays visible
for perhaps 1e4 to a few 1e5 years before it dissolves into the interstellar
medium. A runaway at 1,000-2,500 km/s travels 1-2.5 kpc per Myr, so only the
youngest survivors can still be near theirs. Shen et al. (2018) found that
D6-2's path passes 0.9 deg from the centre of the faint, old remnant
G70.0-21.5 about 90,000 years ago, and left the chance of a coincidence "to
be quantified in future work". We quantify it here.

The geometry. The binary orbited with the Galactic disc before it exploded,
and the remnant's centre keeps that motion. So the site of an explosion that
happened a time t ago, *as it appears today*, is the star's position back then
plus the disc's circular velocity there times t:

    remnant_today(t) = x_star(-t) + v_rot(x_star(-t)) * t

Each posterior sample gives such a curve; seen from the Sun, it is a track
across the sky that starts at the star. A remnant matches if its centre lies
within its own radius of the track, at a time t younger than the look-back
limit.

The chance test. Keeping every remnant's latitude and size but shifting its
longitude at random destroys any physical connection while preserving how
remnants crowd toward the Galactic plane. Repeating that thousands of times
gives the probability that a star's track hits *some* remnant by chance as
closely as the real one does.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import Galactic, SkyCoord
from astropy.table import Table

from config import AUX_DIR
from kinematics import galcen_frame, rotation_velocity

GREEN_PATH = AUX_DIR / "green2025_snrs.tsv"
YR_PER_MYR = 1e6


def load_green() -> Table:
    """Green (2025): 310 remnants, with Galactic coordinates and sizes."""
    rows = []
    for ln in GREEN_PATH.read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        parts = ln.split("\t")
        if not parts[0].strip() or not parts[0].strip()[0].isdigit():
            continue
        l_, b_ = float(parts[0]), float(parts[1])
        maj = float(parts[6]) if parts[6].strip() else np.nan
        mn = float(parts[7]) if parts[7].strip() else maj
        rows.append(dict(name=parts[3].strip(), l=l_, b=b_, maj_arcmin=maj, min_arcmin=mn,
                         type=parts[9].strip(), names=parts[15].strip()))
    t = Table(rows=rows)
    # radius: half the geometric-mean diameter
    t["radius_deg"] = 0.5 * np.sqrt(t["maj_arcmin"] * t["min_arcmin"]) / 60.0
    return t


def apparent_site_track(pos, t_myr, vcirc_fn):
    """Today's apparent position of an explosion site, for explosions at each
    past time along the orbits. pos: (3, n, nt) over t_myr (0 -> past).
    Returns l, b (deg) and heliocentric distance (kpc), each (n, nt)."""
    n, nt = pos.shape[1], pos.shape[2]
    flat = pos.reshape(3, -1)
    vrot = rotation_velocity(flat, vcirc_fn)                    # km/s
    age_myr = np.broadcast_to(-np.asarray(t_myr)[None, :], (n, nt)).reshape(-1)
    site = flat + vrot * age_myr * 1.0227121650537077e-3        # kpc
    g = SkyCoord(x=site[0] * u.kpc, y=site[1] * u.kpc, z=site[2] * u.kpc,
                 frame=galcen_frame()).transform_to(Galactic())
    return (g.l.deg.reshape(n, nt), g.b.deg.reshape(n, nt),
            g.distance.to_value(u.kpc).reshape(n, nt))


def _sep_deg(l1, b1, l2, b2):
    l1, b1, l2, b2 = map(np.radians, (l1, b1, l2, b2))
    s = (np.sin((b2 - b1) / 2) ** 2 + np.cos(b1) * np.cos(b2) * np.sin((l2 - l1) / 2) ** 2)
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(s, 0, 1))))


def closest_approach(l_track, b_track, t_yr, rem_l, rem_b):
    """For one track (nt,) and remnants (m,): minimum separation (deg) and the
    look-back time (yr) at which it happens."""
    sep = _sep_deg(l_track[None, :], b_track[None, :], np.asarray(rem_l)[:, None],
                   np.asarray(rem_b)[:, None])
    k = sep.argmin(axis=1)
    return sep[np.arange(len(k)), k], np.asarray(t_yr)[k]


def match_samples(l, b, dist, t_yr, green: Table, max_age_yr: float):
    """For every remnant, the fraction of posterior tracks that pass inside its
    radius within max_age_yr, and the median closest approach and time.

    l, b, dist: (n, nt) apparent-site tracks; t_yr: (nt,) look-back times."""
    keep = np.asarray(t_yr) <= max_age_yr
    l, b, t_yr = l[:, keep], b[:, keep], np.asarray(t_yr)[keep]
    rows = []
    # cheap pre-filter: remnants within 30 deg of any track point
    lc, bc = np.median(l, axis=0), np.median(b, axis=0)
    near = np.zeros(len(green), bool)
    for k in range(0, len(lc), max(1, len(lc) // 30)):
        near |= _sep_deg(lc[k], bc[k], green["l"], green["b"]) < 30
    for j in np.flatnonzero(near):
        seps = _sep_deg(l, b, green["l"][j], green["b"][j])           # (n, nt)
        kmin = seps.argmin(axis=1)
        smin = seps[np.arange(len(kmin)), kmin]
        inside = smin < green["radius_deg"][j]
        rows.append(dict(name=green["name"][j], radius_deg=float(green["radius_deg"][j]),
                         frac_inside=float(inside.mean()), sep_median=float(np.median(smin)),
                         t_median_yr=float(np.median(t_yr[kmin])),
                         d_star_median=float(np.median(dist[np.arange(len(kmin)), kmin]))))
    rows.sort(key=lambda r: (-r["frac_inside"], r["sep_median"]))
    return rows


def chance_probability(l_med, b_med, green: Table, rng, n_shuffles: int,
                       best_sep_over_radius: float) -> dict:
    """Probability that a track (median, (nt,)) passes some remnant at least as
    deep inside it (separation / radius) as the real best match, when every
    remnant's longitude is shifted at random and latitudes are kept."""
    rad = np.asarray(green["radius_deg"])
    gl, gb = np.asarray(green["l"]), np.asarray(green["b"])
    # only remnants whose latitude band the track crosses can ever match
    bmin, bmax = b_med.min(), b_med.max()
    cand = (gb + rad >= bmin - 1) & (gb - rad <= bmax + 1)
    gl, gb, rad = gl[cand], gb[cand], rad[cand]
    hits = 0
    best = np.empty(n_shuffles)
    for k in range(n_shuffles):
        shifted = (gl + rng.uniform(0, 360, size=gl.size)) % 360
        sep = _sep_deg(l_med[None, :], b_med[None, :], shifted[:, None], gb[:, None]).min(axis=1)
        ratio = (sep / rad).min() if sep.size else np.inf
        best[k] = ratio
        hits += ratio <= best_sep_over_radius
    return dict(p_chance=hits / n_shuffles, n_shuffles=n_shuffles,
                n_remnants_in_band=int(cand.sum()),
                median_best_ratio=float(np.median(best)))
