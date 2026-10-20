"""Published values, used only for comparison and labelled as such in every figure.

Each entry was checked against the paper itself (abstract or text) on 2026-10-04.
"""

from __future__ import annotations

# --- J0313-1806, the most distant quasar known when found (Wang et al. 2021, ApJL 907, L1) ---
WANG21_J0313 = dict(
    z=7.6423, z_err=0.0013,              # [C II] redshift from ALMA
    mbh=1.6e9, mbh_err=0.4e9,            # Mg II, Vestergaard & Osmer 2009 recipe
    fwhm_mgii=3670.0, fwhm_mgii_err=405.0,
    lambda_l3000=2.7e46,                 # erg/s (their Table 1)
    lbol=1.4e47, lbol_lsun=3.6e13,
    edd_ratio=0.67, edd_ratio_err=0.14,
    age_myr=670.0,                       # "~670 million years after the Big Bang"
    seed_range=(1e4, 1e5),               # Eddington growth from z ~ 15-30
    seed_z_start=(15.0, 30.0),
)

# --- Vanden Berk et al. 2001 (AJ 122, 549): composite of 2,204 SDSS quasars ---
VB01 = dict(
    n_spectra=2204,
    alpha_nu_gm=-0.44,        # geometric mean, windows 1350-1365 and 4200-4230 A
    alpha_nu_gm_red=-2.45,    # geometric mean, 6005-6035 and 7160-7180 A
    alpha_nu_median=-0.46,
    alpha_nu_median_red=-1.58,
    alpha_nu_err=0.1,         # dominated by spectrophotometric calibration
    wave_range=(800.0, 8555.0),
)
#: Richards et al. 2006: their DR1 composite, alpha_nu = -0.436 from Ly alpha to H beta.
R06_ALPHA_NU = -0.436

# --- Luminosity function ---
R06 = dict(
    n_quasars=15343, area_deg2=1622.0,
    luminous_cut=-27.6,       # complete within the flux limits at 0.8 < z < 4.8
    peak_statement="the space density of luminous quasars peaks between redshifts 2 and 3",
)
SK12 = dict(n_quasars=57959, area_deg2=6248.0)

# --- Eddington ratios ---
#: Kollmeier et al. 2006 (ApJ 648, 128): 407 AGES AGNs, 0.3 < z < 4, L_bol 1e45-1e47:
#: log-normal with peak L/L_Edd ~ 1/4 and dispersion 0.3 dex.
KOLLMEIER06 = dict(peak=0.25, sigma_dex=0.3, n=407)

# --- C IV masses ---
#: Coatman et al. 2017 (MNRAS 465, 2120): uncorrected C IV masses are too high by a
#: factor of more than five at blueshifts of 3000 km/s.
COATMAN17_BIAS_AT_3000 = 5.0

# --- Black hole growth ---
#: Typical systematic uncertainty of single-epoch virial masses (Shen 2013, BASI 41, 61).
SE_MASS_SYSTEMATIC_DEX = 0.4

# --- Cosmic star formation, for shape comparison (Madau & Dickinson 2014, eq. 15) ---
def md14_sfrd(z):
    """psi(z) in Msun/yr/Mpc^3 (Salpeter IMF)."""
    return 0.015 * (1 + z) ** 2.7 / (1 + ((1 + z) / 2.9) ** 5.6)


# --- Black holes found by JWST (published values; methods differ, see notes) ---
JWST_BLACK_HOLES = [
    dict(name="GN-z11", z=10.603, log_mbh=6.2, log_mbh_err=0.3, edd_ratio=5.0,
         method="UV lines, local virial relations", ref="Maiolino et al. 2024, Nature 627, 59"),
    dict(name="UHZ1", z=10.073, log_mbh=7.5, log_mbh_err=0.5, edd_ratio=1.0,
         method="X-ray luminosity, Eddington-limited accretion assumed (1e7-1e8)",
         ref="Bogdan et al. 2024, Nat. Astron. 8, 126; z: Goulding et al. 2023"),
    dict(name="CEERS_1019", z=8.679, log_mbh=6.95, log_mbh_err=0.37, edd_ratio=1.2,
         method="broad H beta", ref="Larson et al. 2023, ApJL 953, L29"),
    dict(name="Abell2744-QSO1", z=7.0451, log_mbh=7.6, log_mbh_err=0.15, edd_ratio=0.3,
         method="broad H beta, a 'little red dot'", ref="Furtak et al. 2024, Nature 628, 57"),
]
