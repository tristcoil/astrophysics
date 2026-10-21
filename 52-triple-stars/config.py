"""Project configuration: paths, data sources, published system parameters.

Every choice that changes a number quoted in the README lives here, once:
which archives and catalogue versions we read, which TESS sectors, the
published orbital elements that seed the simulations (with their sources),
quality cuts and random seeds.

Conventions used throughout the project
---------------------------------------
* Masses in solar masses, lengths in au (or solar radii where stated), times in
  days for light curves and Julian years for the long-term dynamics.
* Orbital elements follow the usual binary-star conventions: ``omega`` is the
  argument of periastron of the *secondary* (the companion), ``Omega`` the
  position angle of the ascending node, ``i`` the inclination of the orbit to
  the plane of the sky (i < 90 deg: counter-clockwise on the sky).
* The mutual inclination J of two orbits follows
  cos J = cos i1 cos i2 + sin i1 sin i2 cos(Omega1 - Omega2).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RAW = DATA / "raw"              # downloads, never modified
DERIVED = DATA / "derived"      # parsed tables, light curves, eclipse times
CACHE = DATA / "cache"          # expensive results: fits, N-body runs
FIGURES = HERE / "figures"
MANIFEST = DATA / "manifest.json"
RESULTS = DATA / "results.json"

for _d in (DATA, RAW, DERIVED, CACHE, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 52003         # fixed seed: every Monte Carlo draw is reproducible
N_PROCESSES = 8             # worker processes for N-body grids (machine is shared)

# ---------------------------------------------------------------------------
# Physical constants (IAU 2015 nominal values, as astropy uses them)
# ---------------------------------------------------------------------------

GM_SUN = 1.3271244e20            # m^3 s^-2
R_SUN = 6.957e8                  # m
AU = 1.495978707e11              # m
C_LIGHT = 299792458.0            # m s^-1
DAY = 86400.0                    # s
YEAR = 365.25 * DAY              # Julian year, s
R_SUN_AU = R_SUN / AU
#: G in au^3 Msun^-1 day^-2 and au^3 Msun^-1 yr^-2 (from the constants above)
G_AU_DAY = GM_SUN * DAY**2 / AU**3
G_AU_YR = GM_SUN * YEAR**2 / AU**3
#: speed of light in au/day and au/yr (REBOUNDx 'gr' needs c in sim units)
C_AU_DAY = C_LIGHT * DAY / AU
C_AU_YR = C_LIGHT * YEAR / AU
#: light-travel time across 1 au, in days
LIGHT_AU_DAY = AU / C_LIGHT / DAY
JD_J2000 = 2451545.0


def jd_to_year(jd):
    """Julian epoch (decimal year) from a Julian date."""
    return 2000.0 + (np.asarray(jd, float) - JD_J2000) / 365.25


def year_to_jd(year):
    return JD_J2000 + (np.asarray(year, float) - 2000.0) * 365.25


def semimajor_axis_au(m_tot_msun: float, period_days: float) -> float:
    """Kepler's third law with this project's G."""
    return (G_AU_DAY * m_tot_msun * period_days**2 / (4 * np.pi**2)) ** (1 / 3)


# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: Tokovinin's Multiple Star Catalog, latest release (19 June 2026; 6,080
#: systems). Plain http: the site does not serve https.
MSC_VERSION = "20260619"
MSC_URL = f"http://www.ctio.noirlab.edu/~atokovin/stars/newmsc-{MSC_VERSION}.tar.gz"
MSC_README_URL = "http://www.ctio.noirlab.edu/~atokovin/stars/Readme"

#: CDS (VizieR) tables used for the statistics.
CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp"
VIZIER_TABLES = {
    # Tokovinin 2017, ApJ 844, 103: orbit alignment in triple stars
    "J/ApJ/844/103": ("ReadMe", "table1.dat", "table3.dat", "refs.dat"),
    # Tokovinin, Thomas, Sterzik & Udry 2006, A&A 450, 681: tertiary companions
    # to 165 close spectroscopic binaries
    "J/A+A/450/681": ("ReadMe", "table2.dat", "table5.dat", "notes.dat"),
    # Claret 2017, A&A 600, A30: quadratic limb darkening for TESS (ATLAS)
    "J/A+A/600/A30": ("ReadMe", "table25.dat"),
    # Claret 2000, A&A 363, 1081: limb darkening for UBVRI (ATLAS); B band for the plates
    "J/A+A/363/1081": ("ReadMe", "atlasco.dat"),
}

#: DASCH DR7 (Harvard plates, 1885-1992) public web API. No key needed; the
#: light curves are calibrated against APASS (B band for blue plates).
DASCH_API = "https://api.starglass.cfa.harvard.edu/public"
DASCH_REFCAT = "apass"
#: daschlab's "standard rejections" (apply_standard_rejections, v1.0.0):
#: AFLAGS bits HIGH_BACKGROUND, LARGE_ISO_RMS, LARGE_LOCAL_SMOOTH_RMS,
#: CLOSE_TO_LIMITING, BIN_DRAD_UNKNOWN.
DASCH_BAD_AFLAGS = (1 << 6) | (1 << 11) | (1 << 12) | (1 << 13) | (1 << 15)
#: We also drop BAD_PLATE_QUALITY (1 << 7), UNCERTAIN_DATE (1 << 9: an eclipse
#: needs a good time stamp) and MULT_EXP_BLEND (1 << 10).
DASCH_EXTRA_BAD_AFLAGS = (1 << 7) | (1 << 9) | (1 << 10)
#: largest acceptable time uncertainty of a plate, as a fraction of the
#: eclipse duration of each system (set per system below)
DASCH_TIME_ACCURACY_FRACTION = 0.1

#: MAST services
MAST_TESSCUT = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut"
MAST_DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file?uri="
TESSCUT_SIZE = 15               # pixels (21" each)

#: Gaia DR3 (ESA archive) for Proxima Centauri
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
GAIA_TAP_MIRROR = "https://gaia.ari.uni-heidelberg.de/tap/sync"
PROXIMA_GAIA_DR3 = 5853498713190525696

# ---------------------------------------------------------------------------
# The systems. Every number carries its source; these seed the simulations
# and define the eclipse models. Masses Msun, radii Rsun, periods days.
# ---------------------------------------------------------------------------

SYSTEMS = {
    # ------------------------------------------------------------------
    "SS Lac": dict(
        wds="22047+4626", tic=468682885, ra=331.173131, dec=46.427252,
        note="Eclipsed ~1890-1950, then stopped (NGC 7209 member).",
        # inner pair: Torres & Stefanik 2000 (AJ 119, 1914) via SB9 2456 / MSC;
        # sidereal period from TS00 as quoted by Torres 2001 (AJ 121, 2227)
        P_in=14.4161471, P_in_err=0.0000089, e_in=0.136, omega_in=178.28,
        T_peri_in=2450381.0742, K1=74.69, K2=76.77,
        # outer orbit of the binary around the barycentre (TS00, SB9 2458)
        P_out=678.8, e_out=0.159, omega_out=256.70, T_peri_out=2450324.0, K_AB=5.41,
        # Torres 2001: eclipses possible for 81.6 <= i <= 98.4 deg (TS00 light
        # curve), so (r1 + r2)(1 - e^2)^-1 ~ cos(81.6 deg)
        i_eclipse_limit=81.6,
        # Torres 2001 nodal-model fit (P_node fixed at 600 yr)
        pub=dict(P_node=600.0, I=84.0, I_err=10.0, eps=15.0, eps_err=3.0, t90=1905.8, t90_err=1.5,
                 didt=-0.1441, didt_err=0.0081, t90_lin=1905.6, t90_lin_err=1.4,
                 J_range=(24.0, 30.0), end_primary=1945, end_secondary=1950,
                 next_eclipses="shortly after 2200",
                 ref="Torres 2001, AJ 121, 2227 (astro-ph/0012542); TS00"),
        eclipse_hw_days=0.30,      # ~half-width of an eclipse (0.02 in phase)
    ),
    # ------------------------------------------------------------------
    "HS Hya": dict(
        wds="10246-1906", tic=434479378, ra=156.153199, dec=-19.092488,
        note="Eclipses appeared ~1920, deepest ~1960, gone by 2021.",
        # Torres et al. 1997 (A&A 318, 117) masses/radii/temperatures as quoted
        # by Davenport et al. 2021 (Table 1); period Popper 1971; epoch
        # Gyldenkerne et al. 1975
        P_in=1.568024, e_in=0.0, T0_in=2441374.5934,
        m1=1.255, m2=1.219, R1=1.2747, R2=1.2161, Teff1=6500.0, Teff2=6400.0,
        # outer orbit: Torres et al. 1997; MSC elements from Vokrouhlicky & Zasche 2022
        P_out=190.53, e_out=0.246, omega_out=111.2, T_peri_out=2448047.2, K_AB=9.02,
        m3_guess=0.56,
        pub=dict(m1=1.31, m1_err=0.03, m2=1.27, m2_err=0.03, m3=0.56, m3_lo=0.47, m3_hi=0.68,
                 J_prograde=(50.0, 65.0), J_retrograde=(120.0, 150.0), P_node=700.0,
                 P_node_err=200.0, next_eclipses=2200,
                 ref="Vokrouhlicky & Zasche 2022, AJ 163, 94",
                 # Davenport et al. 2021 (AJ 162, 189): TESS S9 primary depth
                 dav_depth_s9_mag=0.00173, dav_depth_s9_err=0.00007, dav_return=2195,
                 dav_tumble=1194.0, dav_end=2021.1),
        eclipse_hw_days=0.12,
    ),
    # ------------------------------------------------------------------
    "V907 Sco": dict(
        wds="17569-3445", tic=261030551, ra=269.231509, dec=-34.750338,
        note="Eclipsing before ~1920, in 1963-1986, and again since 2019.",
        # Zasche, Vokrouhlicky, Barlow & Masek 2023 (AJ 165, 81), Table 1
        P_in=3.7762049, P_in_err=0.0000122, e_in=0.0080, omega_in=122.6,
        T0_in=2459362.3272, a_in_rsun=17.79,
        m1=2.74, m2=2.56, R1=2.68, R2=2.53, Teff1=10000.0, Teff2=10004.0,
        L3_tess=0.059,
        # outer orbit (Zasche et al. 2023): P2 = 142.01 d, e2 = 0.47, m3 = 1.06
        P_out=142.01, e_out=0.47, m3_guess=1.06,
        pub=dict(J=26.2, J_lo=24.0, J_hi=28.8, i_tot=71.0, P_node=63.5, P_node_lo=60.9,
                 P_node_hi=66.8, eclipsing_until=2049, m3=1.06, m3_lo=0.96, m3_hi=1.17,
                 ref="Zasche et al. 2023, AJ 165, 81 (arXiv:2302.04613)",
                 # their Table 2 (mid HJD - 2400000, i1, sigma)
                 inclinations=[(15207.0, 86.0, 7.3, "DASCH"), (16835.5, 86.6, 2.9, "DASCH"),
                               (18312.0, 88.0, 1.8, "DASCH"), (19787.5, 87.3, 1.6, "DASCH"),
                               (20910.0, 83.7, 3.4, "DASCH"), (21793.0, 81.2, 3.7, "DASCH"),
                               (38420.5, 86.23, 1.85, "APPLAUSE"),
                               (58675.0, 73.37, 0.25, "TESS S13"),
                               (59365.3, 76.57, 0.02, "TESS S39a"),
                               (59383.9, 76.61, 0.02, "TESS S39b"),
                               (59673.8, 77.78, 0.10, "FRAM"), (59771.0, 78.11, 0.08, "FRAM"),
                               (59794.7, 78.45, 0.04, "FRAM"), (59812.0, 78.55, 0.04, "FRAM")],
                 i_limit=73.3),
        eclipse_hw_days=0.13,
    ),
    # ------------------------------------------------------------------
    "TIC 290061484": dict(
        wds="21044+5114", tic=290061484, ra=316.109369, dec=51.225478,
        gaia_dr3=2169382208774963072,
        note="Triply eclipsing; record-short 24.5-d outer orbit; (2+1)+1 quadruple.",
        # Kostov et al. 2024, ApJ 974, 25 (arXiv:2410.01711), Table 3:
        # osculating elements at t0 = BJD 2458711.0
        t0_bjd=2458711.0,
        Aa_Ab=dict(P=1.79364, a_rsun=14.59, e=0.00250, omega=151.0, i=85.50, Omega=0.0,
                   T_inf=2458713.0340, tau=2458712.442),
        A_B=dict(P=24.6347, a_rsun=98.18, e=0.2011, omega=92.2, i=94.43, Omega=1.70,
                 T_inf=2458728.733, tau=2458704.228),
        AB_C=dict(P=3205.0, a_rsun=2734.0, e=0.44, omega=155.0, i=82.0, Omega=-1.0,
                  tau=2458931.0),
        masses=dict(Aa=6.853, Ab=6.106, B=7.903, C=6.009),
        radii=dict(Aa=3.077, Ab=2.803, B=3.499, C=2.779),
        teff=dict(Aa=22154.0, Ab=20887.0, B=23703.0, C=20690.0),
        flux_tess=dict(Aa=0.2478, Ab=0.1866, B=0.3601, C=0.1804),
        pub=dict(P_out_mean=24.498, P_out_mean_err=0.005, P_out_linear=24.489, P_out_linear_err=0.004,
                 i_mut=9.17, i_mut_err=0.4, P_node_yr=2.40, P_apse_dyn_yr_in=1.34,
                 ltte_amp_triple_d=0.0009, distance_pc=1519, P_in_over_P_out=13.7,
                 tertiary_events=[(8728.46, "primary"), (8753.25, "primary"), (9831.33, "primary"),
                                  (9840.48, "secondary"), (10344.80, "primary"),
                                  (10353.98, "secondary"), (10369.38, "primary"),
                                  (10378.60, "secondary"), (10393.89, "primary")],
                 ref="Kostov et al. 2024, ApJ 974, 25"),
        sectors=(15, 16, 55, 56, 75, 76, 82, 83),
    ),
    # ------------------------------------------------------------------
    "Algol": dict(
        wds="03082+4057", ra=47.042209, dec=40.955639,
        note="Inner and outer orbits perpendicular (CHARA).",
        # Baron et al. 2012, ApJ 752, 20 (Tables 4-6)
        m_A=3.17, m_B=0.70, m_C=1.76, R_A=2.73, R_B=3.48, R_C=1.73,
        inner=dict(P=2.867328, a_mas=2.15, i=98.70, Omega=43.43, e=0.0, T0=2441771.353),
        outer=dict(P=680.168, a_mas=93.43, i=83.66, Omega=132.66, e=0.227, omega=310.02,
                   T0=2446927.22),
        parallax_mas=34.7,
        pub=dict(J=90.20, J_err=0.32, ref="Baron et al. 2012, ApJ 752, 20"),
    ),
    # ------------------------------------------------------------------
    "alpha Cen": dict(
        wds="14396-6050",
        note="The nearest triple: alpha Cen A + B (80 yr) with Proxima (~550 kyr).",
        # Kervella, Thevenin & Lovis 2017 (A&A 598, L7), Tables 1-3: inputs at
        # the Hipparcos epoch J1991.25; RVs corrected for gravitational redshift
        # and convective blueshift
        k17=dict(
            epoch=1991.25,
            AB=dict(ra="14:39:40.2068", ra_err_mas=25.0, dec="-60:50:13.673", dec_err_mas=19.0,
                    plx=747.17, plx_err=0.61, pmra=-3619.9, pmra_err=3.9, pmdec=693.8,
                    pmdec_err=3.9, rv=-22.332, rv_err=0.005),
            Prox=dict(ra="14:29:47.7474", ra_err_mas=1.3, dec="-62:40:52.868", dec_err_mas=1.5,
                      plx=768.77, plx_err=0.37, pmra=-3773.8, pmra_err=0.4, pmdec=770.5,
                      pmdec_err=2.0, rv=-22.204, rv_err=0.032),
            m_A=1.1055, m_A_err=0.0039, m_B=0.9373, m_B_err=0.0033, m_C=0.1221, m_C_err=0.0022,
            out=dict(a_kau=8.7, a_lo=8.3, a_hi=9.4, e=0.50, e_lo=0.41, e_hi=0.58, P_kyr=547,
                     P_lo=507, P_hi=613, i=107.6, i_lo=105.6, i_hi=109.4, Omega=126.0, Omega_err=5.0,
                     omega=72.3, omega_lo=65.7, omega_hi=81.0, T0_kyr=283, T0_lo=242, T0_hi=342,
                     peri_kau=4.3, peri_lo=3.4, peri_hi=5.4, apo_kau=13.0, apo_lo=12.9, apo_hi=13.3,
                     v_rel=273.0, v_rel_err=49.0, v_esc=545.0, v_esc_err=11.0, d_au=12947.0,
                     d_au_err=260.0, mutual=30.0),
            ref="Kervella, Thevenin & Lovis 2017, A&A 598, L7"),
        # Akeson et al. 2021 (AJ 162, 14), Tables 8-9: ALMA + Hipparcos, J2019.5
        a21=dict(epoch=2019.5, ra_deg=219.85892215, ra_err_mas=1.50, dec_deg=-60.83163195,
                 dec_err_mas=1.17, pmra=-3639.95, pmra_err=0.42, pmdec=700.40, pmdec_err=0.17,
                 plx=750.81, plx_err=0.38, V0=-22.3796, V0_err=0.0020, m_A=1.0788, m_A_err=0.0029,
                 m_B=0.9092, m_B_err=0.0025,
                 AB_orbit=dict(P_yr=79.762, a_arcsec=17.4930, e=0.51947, i=79.2430, Omega=205.073,
                               omega=231.519, T0=1955.564),
                 ref="Akeson et al. 2021, AJ 162, 14"),
        # gravitational-redshift difference alpha Cen A - Sun applied by K17 to
        # the barycentric RV of AB (m/s)
        dv_grs_A=-61.4, dv_grs_A_err=2.7,
    ),
    # ------------------------------------------------------------------
    "HD 181068": dict(
        wds="19171+4116", kic=5952403, ra=289.287411, dec=41.264810,
        note="Red giant eclipsed by, and eclipsing, a pair of red dwarfs (Kepler).",
        P_in=0.9057, P_out=45.4711, ref="Derekas et al. 2011, Science 332, 216; Borkovits et al. 2013"),
    "KOI-126": dict(
        wds="19499+4107", kic=5897826, ra=297.475862, dec=41.114260,
        note="Two M dwarfs (1.77 d) transiting a 1.35-Msun star every 33.9 d (Kepler).",
        P_in=1.76713, P_out=33.9214, masses=(1.347, 0.2413, 0.2127),
        ref="Carter et al. 2011, Science 331, 562"),
}

#: TESS sectors fetched as full-frame-image cutouts (own photometry) or SPOC
#: 2-minute light curves.
TESS_PLAN = {
    "TIC 290061484": dict(kind="tesscut", sectors=(15, 16, 55, 56, 75, 76, 82, 83)),
    "V907 Sco": dict(kind="tesscut", sectors=(13, 39, 66, 91, 92, 93)),
    "SS Lac": dict(kind="tesscut", sectors=(16, 56, 76, 83)),
    "HS Hya": dict(kind="spoc", sectors=(9, 35, 62, 89, 99, 100)),
}
#: Kepler long-cadence quarters for the triply-eclipsing gallery
KEPLER_PLAN = {
    "HD 181068": dict(kic=5952403, quarters=(1, 2, 3, 4, 5, 6)),
    "KOI-126": dict(kic=5897826, quarters=(1, 2, 3, 4, 5, 6)),
}

# ---------------------------------------------------------------------------
# Stability criterion (Mardling & Aarseth 2001, MNRAS 321, 398)
# ---------------------------------------------------------------------------

MA01_C = 2.8


def ma01_ratio(q_out, e_out, i_mut_deg=0.0):
    """Critical a_out/a_in of Mardling & Aarseth (2001):

        a_out/a_in > 2.8 (1+q_out)^(2/5) (1+e_out)^(2/5) (1-e_out)^(-6/5) (1 - 0.3 i/180 deg)

    with q_out = m3/(m1+m2) and i the mutual inclination in degrees (their
    empirical inclination factor, 1 for prograde coplanar, 0.7 for retrograde).
    """
    q_out = np.asarray(q_out, float)
    e_out = np.asarray(e_out, float)
    return (MA01_C * (1 + q_out) ** 0.4 * (1 + e_out) ** 0.4 * (1 - e_out) ** (-1.2)
            * (1 - 0.3 * np.asarray(i_mut_deg, float) / 180.0))


def ma01_period_ratio(m_in, m3, e_out, i_mut_deg=0.0):
    """The same limit as a ratio of periods: P_out/P_in = (a_out/a_in)^1.5
    sqrt(m_in / (m_in + m3))."""
    m_in = np.asarray(m_in, float)
    m3 = np.asarray(m3, float)
    x = ma01_ratio(m3 / m_in, e_out, i_mut_deg)
    return x**1.5 * np.sqrt(m_in / (m_in + m3))
