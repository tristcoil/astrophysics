"""Load the cached datasets into tidy tables, with every cut documented.

The CHIME catalogue has one row per *sub-burst*: a burst with three components
gets three rows that share its position, DM and fluence. Bursts in turn belong to
*sources* -- most sources have been seen once, but 83 are repeaters with up to
hundreds of bursts. Counting the wrong unit is the easiest mistake to make with
this catalogue, so the three levels are kept distinct here:

    rows (5045)  ->  bursts: sub_num == 0 (4539)  ->  sources: first burst each (3641)

Population statistics use one entry per source, taken from its first burst, and
drop events CHIME flags for exclusion -- the convention of the catalogue paper.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import Galactic, SkyCoord
from astropy.io import fits
from astropy_healpix import HEALPix

from config import (
    ATNF_PATH,
    BAND_EDGE_TOL_MHZ,
    BAND_HIGH_MHZ,
    BAND_LOW_MHZ,
    CHIME_CATALOG_PATH,
    DM_HALO,
    FRB_REPO_DIR,
    MORPHOLOGY_FLAGS,
    P_OX_MIN,
    POPULATION_FLAGS,
    Z_MIN,
)

# ---------------------------------------------------------------------------
# CHIME/FRB Catalog 2
# ---------------------------------------------------------------------------


def _require(path):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")
    return path


def load_chime_rows() -> pd.DataFrame:
    """The raw catalogue: one row per sub-burst."""
    return pd.read_csv(_require(CHIME_CATALOG_PATH))


def chime_bursts(rows: pd.DataFrame) -> pd.DataFrame:
    """One row per burst (its first sub-burst), with derived columns.

    Derived quantities, and why:

    * ``source`` -- the repeater's name for repeat bursts, else the burst's own
      TNS name. Groups bursts by the object that emitted them.
    * ``dm_mw_ne2001`` / ``dm_mw_ymw16`` -- the Milky Way disk's DM along the
      sightline under each Galactic electron model, recovered as DM minus the
      catalogue's ``dm_exc_*``. We checked that ``dm_exc_*`` subtracts the disk
      model *only*, no halo term: at |b| > 60 deg the implied NE2001 value has
      median 25 pc cm^-3, exactly the model's polar column.
    * ``bandwidth_mhz`` -- full width at tenth maximum of the fitted spectrum,
      clipped to CHIME's 400-800 MHz. ``band_filling`` marks bursts whose
      emission reaches both band edges: their bandwidth is a lower limit.
    * ``width_ms`` -- the intrinsic Gaussian sigma of the first component, after
      the fit has removed scattering and dispersion smearing.
    """
    n_comp = rows.groupby("tns_name")["sub_num"].max().rename("n_components") + 1
    b = rows[rows["sub_num"] == 0].merge(n_comp, left_on="tns_name", right_index=True)
    b = b.reset_index(drop=True)

    b["is_repeater"] = b["repeater_name"].notna()
    b["source"] = np.where(b["is_repeater"], b["repeater_name"], b["tns_name"])

    b["dm_mw_ne2001"] = b["dm_fitb"] - b["dm_exc_ne2001"]
    b["dm_mw_ymw16"] = b["dm_fitb"] - b["dm_exc_ymw16"]

    b["bandwidth_mhz"] = b["high_freq"] - b["low_freq"]
    b["band_filling"] = (
        (b["low_freq"] <= BAND_LOW_MHZ + BAND_EDGE_TOL_MHZ)
        & (b["high_freq"] >= BAND_HIGH_MHZ - BAND_EDGE_TOL_MHZ)
    )
    b["width_ms"] = b["width_fitb"] * 1e3
    b["scat_ms_400"] = b["scat_time"] * 1e3  # scattering time referenced to 400.2 MHz
    b["exposure_h"] = b["exp_up"].fillna(0) + b["exp_low"].fillna(0)

    b["passes_population"] = (b[list(POPULATION_FLAGS)] == 0).all(axis=1)
    b["has_fit"] = b["dm_fitb"].notna() & b["width_fitb"].notna()
    b["passes_morphology"] = (b[list(MORPHOLOGY_FLAGS)] == 0).all(axis=1) & b["has_fit"]
    return b


def chime_sources(bursts: pd.DataFrame) -> pd.DataFrame:
    """One row per source: its chronologically first burst, plus a burst count."""
    b = bursts.sort_values("mjd_400")
    first = b.groupby("source", sort=False).head(1).copy()
    counts = b.groupby("source").size().rename("n_bursts")
    first = first.merge(counts, left_on="source", right_index=True)
    return first.sort_values("mjd_400").reset_index(drop=True)


def population_sources(sources: pd.DataFrame) -> pd.DataFrame:
    """Sources whose first burst passes CHIME's population-study flags."""
    return sources[sources["passes_population"]].reset_index(drop=True)


def morphology_samples(bursts: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """The three burst samples compared in the repeater figure.

    * ``one_off``: one burst per never-seen-to-repeat source.
    * ``repeater_first``: for each repeater, its earliest burst with a good
      morphology fit -- one entry per source, so a prolific repeater (FRB
      20220912A alone has 373 bursts) cannot swamp the comparison.
    * ``repeater_all``: every good repeat burst, for context only.

    The same flag cuts apply to all three. Unlike the population sample, a
    repeater whose very first burst is flagged is not dropped: we take its first
    *usable* burst, since the question here is burst shape, not detection rate.
    """
    good = bursts[bursts["passes_morphology"]]
    one_off = good[~good["is_repeater"]]
    rep = good[good["is_repeater"]].sort_values("mjd_400")
    return {
        "one_off": one_off.reset_index(drop=True),
        "repeater_first": rep.groupby("source", sort=False).head(1).reset_index(drop=True),
        "repeater_all": rep.reset_index(drop=True),
    }


def chime_cut_flow(rows, bursts, sources, population) -> list[tuple[str, int]]:
    """Counts at each stage of the selection, for the README and the log."""
    return [
        ("catalogue rows (one per sub-burst)", len(rows)),
        ("bursts", len(bursts)),
        ("  from repeaters", int(bursts["is_repeater"].sum())),
        ("sources", len(sources)),
        ("  repeating sources", int(sources["is_repeater"].sum())),
        ("sources passing population flags", len(population)),
        ("  ... with a fitted DM", int(population["dm_fitb"].notna().sum())),
    ]


# ---------------------------------------------------------------------------
# NE2001: the Milky Way's electrons, as an all-sky map
# ---------------------------------------------------------------------------


class Ne2001Map:
    """DM through the whole Galactic disk (NE2001) in any direction.

    A HEALPix map at NSIDE 128 (0.46 deg pixels) from the FRBs/FRB repository.
    Its FITS header omits the coordinate frame; ``analysis.validate_ne2001_map``
    shows the values match CHIME's own per-burst NE2001 column in Galactic
    coordinates (median ratio 1.02) and not in equatorial ones, so we read it as
    Galactic. Values are interpolated bilinearly between pixel centres.
    """

    def __init__(self, path=None):
        path = _require(path or FRB_REPO_DIR / "ne2001_dm_healpix_map_128.fits")
        with fits.open(path) as hdul:
            hdr = hdul[1].header
            self.values = np.asarray(hdul[1].data["T"], dtype=float).ravel()
        self.nside = int(hdr["NSIDE"])
        order = hdr.get("ORDERING", "RING").strip().lower()
        self.hp = HEALPix(nside=self.nside, order=order, frame=Galactic())
        if self.values.size != self.hp.npix:
            raise ValueError("NE2001 map size does not match its NSIDE")

    def dm(self, l_deg, b_deg):
        """NE2001 DM (pc cm^-3) through the Galaxy toward Galactic (l, b)."""
        l = np.asarray(l_deg, dtype=float) * u.deg
        b = np.asarray(b_deg, dtype=float) * u.deg
        return self.hp.interpolate_bilinear_lonlat(l, b, self.values)

    def dm_radec(self, ra_deg, dec_deg):
        c = SkyCoord(ra=np.asarray(ra_deg) * u.deg, dec=np.asarray(dec_deg) * u.deg).galactic
        return self.dm(c.l.deg, c.b.deg)


# ---------------------------------------------------------------------------
# ATNF pulsars
# ---------------------------------------------------------------------------


def load_pulsars() -> pd.DataFrame:
    """Pulsars with position and DM, labelled by where they live.

    ``magellanic`` pulsars sit in the Large or Small Magellanic Cloud -- other
    galaxies, 50-60 kpc away, whose own gas adds DM. ``globular`` ones are in
    Milky Way globular clusters.
    """
    df = pd.read_csv(_require(ATNF_PATH), comment="#", na_values="*", dtype={"assoc": str, "type": str})
    df["assoc"] = df["assoc"].fillna("")
    df["type"] = df["type"].fillna("")
    df["magellanic"] = df["assoc"].str.contains(r"EXGAL:(?:LMC|SMC)", regex=True)
    df["globular"] = df["assoc"].str.contains("GC:", regex=False)
    return df


# ---------------------------------------------------------------------------
# Localized FRBs with host-galaxy redshifts
# ---------------------------------------------------------------------------

#: Name differences between the repository's two tables. FRB 20180924B (the
#: first FRB localized by ASKAP) is labelled with suffix "A" in FRBs_base.csv;
#: the host table and the Transient Name Server call it "B".
_NAME_ALIASES = {"FRB20180924A": "FRB20180924B"}


def _survey(telescope: str) -> str:
    return {
        "ASKAP": "ASKAP (CRAFT)",
        "DSA": "DSA-110",
        "CHIME": "CHIME Outriggers",
        "MeerKAT": "MeerKAT",
        "VLA": "Other",
    }.get(telescope, "Other")


def load_localized_frbs(ne2001: Ne2001Map | None = None, dm_halo: float = DM_HALO) -> tuple[pd.DataFrame, list]:
    """Localized FRBs with secure host redshifts, and the cut flow that made them.

    Joins the repository's burst table (DM, position, redshift) with its host
    table (PATH association probability P(O|x), host redshift). Where both give
    a redshift we take the host table's: it is the galaxy's measured redshift,
    and the two disagree for only two bursts (the README lists them).

    The extragalactic DM is what remains after removing the Milky Way:

        DM_EG = DM_obs - DM_ISM(NE2001) - DM_halo
    """
    base = pd.read_csv(_require(FRB_REPO_DIR / "FRBs_base.csv"))
    hosts = pd.read_csv(_require(FRB_REPO_DIR / "public_hosts.csv"), dtype={"FRB": str})

    base["name"] = base["Name"].replace(_NAME_ALIASES)
    hosts["name"] = "FRB" + hosts["FRB"].str.strip()
    df = base.merge(hosts[["name", "P_Ox", "z"]].rename(columns={"z": "z_host"}), on="name", how="left")

    df["p_ox"] = df["P_Ox"].combine_first(df["P(O|x)"])
    df["z_base"] = df["z"]
    df["z"] = df["z_host"].combine_first(df["z_base"])
    df["z_disagree"] = (df["z_host"] - df["z_base"]).abs() > 1e-3
    df["survey"] = df["telescope"].map(_survey)
    df["repeater"] = df["repeater"].astype(bool)

    ne2001 = ne2001 or Ne2001Map()
    df["dm_ism"] = ne2001.dm_radec(df["ra"].to_numpy(), df["dec"].to_numpy())
    df["dm_halo"] = dm_halo
    df["dm_eg"] = df["DM"] - df["dm_ism"] - dm_halo

    flow = [("localized FRBs in the compilation", len(df))]
    df = df[df["z"].notna()]
    flow.append(("  with a host redshift", len(df)))
    df = df[df["p_ox"] >= P_OX_MIN]
    flow.append((f"  secure host, P(O|x) >= {P_OX_MIN}", len(df)))
    df = df[df["z"] >= Z_MIN]
    flow.append((f"  z >= {Z_MIN} (beyond the local volume)", len(df)))
    df = df[df["dm_eg"] > 0]
    flow.append(("  positive extragalactic DM", len(df)))

    cols = ["name", "ra", "dec", "DM", "DM_err", "z", "z_base", "z_host", "z_disagree", "p_ox",
            "repeater", "telescope", "survey", "refs", "dm_ism", "dm_halo", "dm_eg"]
    return df[cols].sort_values("z").reset_index(drop=True), flow
