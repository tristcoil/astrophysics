"""Published values this project is graded against.

Collected here, apart from the analysis, so that it is clear what was measured
and what was looked up. Each value was read from the paper's own text (arXiv
source or PDF), not quoted from memory. None of them feeds into a
measurement; they are consulted only afterwards, for comparison.
"""

# ---------------------------------------------------------------------------
# Cosmic star-formation history: Madau & Dickinson (2014, ARA&A 52, 415)
# ---------------------------------------------------------------------------

#: Their equation 15, Salpeter IMF, M_sun / yr / Mpc^3:
#: psi(z) = A (1+z)^B / (1 + ((1+z)/C)^D)
MD14_FIT = {"A": 0.015, "B": 2.7, "C": 2.9, "D": 5.6}

MD14_CLAIMS = {
    # Abstract: "peaked approximately 3.5 Gyr after the Big Bang, at z ~ 1.9"
    "z_peak": 1.9,
    "t_peak_gyr": 3.5,
    # "declined exponentially at later times, with an e-folding timescale of
    # 3.9 Gyr" (section 5: a fit to the data at z < 1)
    "efold_gyr": 3.9,
    # "Half of the stellar mass observed today was formed before z = 1.3.
    # About 25% formed before the peak ..., and another 25% formed after z = 0.7.
    # Less than ~1% of today's stars formed during the epoch of reionization."
    "z_half_mass": 1.3,
    "frac_before_peak": 0.25,
    "z_last_quarter": 0.7,
    "frac_reionization_max": 0.01,
    # Section 6: "rho_* = 5.8 x 10^8 M_sun Mpc^-3 predicted by our model SFH",
    # vs Gallazzi et al. (2008): (6.0 +/- 1.0) x 10^8 (Salpeter).
    "rho_star_z0": 5.8e8,
    "rho_star_gallazzi": (6.0e8, 1.0e8),
    # "mass-weighted mean stellar age ... <t_age> = 8.3 Gyr"
    "mean_stellar_age_gyr": 8.3,
    # "a stellar baryon fraction of only 9% (5% for a Chabrier IMF)"
    "stellar_baryon_fraction": 0.09,
    "stellar_baryon_fraction_chabrier": 0.05,
    # "psi(z) ~ (1+z)^-2.9 at 3 < z < 8"
    "high_z_slope": -2.9,
    # "The comoving SFRD at redshift 7 was approximately the same as that
    # measured locally."
}

#: Cosmology Madau & Dickinson assume: (Omega_M, Omega_L, Omega_b, h).
MD14_COSMOLOGY = {"Om0": 0.3, "Ode0": 0.7, "Ob0": 0.045, "h": 0.7}


# ---------------------------------------------------------------------------
# The local galaxy population
# ---------------------------------------------------------------------------

#: Renzini & Peng (2015, ApJL 801, L29): the main sequence as the ridge line of
#: the star-forming peak in the 3D SFR-mass-number plot, SDSS DR7 MPA-JHU,
#: 0.02 < z < 0.085, AGN excluded, 1/V_max weighted, 0.2 x 0.2 dex bins.
RENZINI_PENG_2015 = {
    "slope": (0.76, 0.01),
    "intercept": (-7.64, 0.02),  # log SFR at log M = 0
    # "The sigma of its best fit Gaussian is 0.3 dex" (at log M = 10.5)
    "sigma_dex": 0.3,
    # the number x SFR ridge sits ~0.2 dex above the number ridge
    "number_times_sfr_offset_dex": 0.2,
    "n_galaxies": 240_000,
}

#: Peng et al. (2010, ApJ 721, 193), Table 2, SDSS DR7 0.02 < z < 0.085:
#: eps_rho = 1 - exp(-(rho/p1)^p2), eps_m = 1 - exp(-(m/p3)^p4), with rho in
#: units of 1 + delta_5 and m in M_sun; red fraction from (U-B) colour,
#: Chabrier masses. f_red = eps_rho + eps_m - eps_rho eps_m (their eq. 6).
PENG_2010 = {
    "log_p1": (1.84, 0.01),
    "p2": (0.60, 0.01),
    "log_p3": (10.56, 0.01),
    "p4": (0.80, 0.01),
    # main sequence: log sSFR = -10.0 - 0.1 (log m - 10.0) for blue galaxies
    "ms_ssfr_at_1e10": -10.0,
    "ms_ssfr_slope": -0.1,
    # "deviations from this simple separable behavior ... no more than
    # +/-0.2 dex in either variable"
    "separability_tolerance_dex": 0.2,
}

#: Brinchmann et al. (2004, MNRAS 351, 1151), eq. 5: total SFR density at
#: z = 0.1 for a Kroupa IMF, in 1e-2 h_70 M_sun / yr / Mpc^3:
#: 1.915 (+0.02/-0.01 random, +0.14/-0.42 systematic).
BRINCHMANN_2004_SFRD = {
    "value": 1.915e-2, "rand": (0.02e-2, 0.01e-2), "sys": (0.14e-2, 0.42e-2),
    "imf": "kroupa", "z": 0.1,
}

#: Salim et al. (2016, ApJS 227, 2), section 7.2: GSWLC vs MPA-JHU (DR7) total
#: specific SFRs for BPT star-forming galaxies agree with scatter 0.18 dex;
#: AGN hosts and line-less galaxies show systematic offsets "in some cases
#: reaching ~1 dex". GSWLC masses exceed MPA-JHU's by 0.03-0.13 dex depending
#: on sSFR, scatter ~0.07 dex.
SALIM_2016 = {
    "ssfr_scatter_sf_dex": 0.18,
    "mass_offset_range_dex": (0.03, 0.13),
    "mass_scatter_dex": 0.07,
}

#: Stellar mass density of the local universe, Salpeter IMF, from Madau &
#: Dickinson's Table 2: Li & White (2009), log rho_* = 8.59 +/- 0.01 at z = 0.07
#: (SDSS); Moustakas et al. (2013), 8.59 +/- 0.05 at 0 < z < 0.2.
LOCAL_STELLAR_MASS_DENSITY = {
    "li_white_2009": (8.59, 0.01),
    "moustakas_2013": (8.59, 0.05),
}

#: Strauss et al. (2002): fibre collisions (55") cost ~6% of Main Galaxy Sample
#: targets; Peng et al. (2010) quote ~10% of targets missed in total.
MGS_FIBRE_LOSS = 0.06
