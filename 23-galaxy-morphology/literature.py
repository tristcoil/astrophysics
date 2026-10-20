"""Published numbers, used only for grading -- never as inputs.

Each entry says which paper, which sample and which definition it refers to,
because a fraction means nothing without the cut that produced it. All values
were read from the papers themselves (arXiv versions).
"""

# Willett et al. 2013, MNRAS 435, 2835 (arXiv:1308.3496): the GZ2 data release.
WILLETT13 = {
    "n_main_specz": 243_500,            # Table 5
    "n_classifications": 16_340_298,    # section 2.4
    "n_volunteers": 83_943,
    "median_classifications": 44,
    # Comparison with Nair & Abraham (2010), section 5.2.
    "na10_overlap": 12_480,
    "na10_oblique_overlap": 5_526,      # p_features > 0.430, p_notedgeon >= 0.715, N >= 20
    "na10_ellipticals_median_psmooth": 0.86,
    "na10_ellipticals_median_psmooth_sd": 0.07,
    "na10_ellipticals_max_pfeatures": 0.741,
    "na10_late_median_pfeatures": 0.88,  # T-types Sa or later
    "smooth08_late_fraction": 0.04,      # p_smooth > 0.8 cut: late types (Sa or later)
    "smooth08_late_or_s0_fraction": 0.20,
    "na10_barred_below_05": 0.39,        # NA10 barred galaxies with GZ2 p_bar < 0.5
    "na10_unbarred_above_05": 0.06,      # NA10 unbarred galaxies with p_bar > 0.5
    "strong_bar_above_08": 0.78, "intermediate_bar_above_08": 0.40,
    "strong_bar_above_05": 0.94, "intermediate_bar_above_05": 0.80,
    "weak_bar_above_08": 0.09, "weak_bar_above_05": 0.32,
}

# Hart et al. 2016, MNRAS 461, 3663 (arXiv:1607.01019): the new debiasing and
# spiral arm number. Luminosity-limited sample 0.03 < z < 0.085, M_r <= -21
# (kcorrect, H0 = 70); mass-limited adds log M* >= 10.6 (Baldry et al. 2006
# masses). Spirals: p_features x p_notedgeon x p_spiral > 0.5 and >= 5 definite
# arm-number votes; arm number = modal debiased fraction. Table 2.
HART16 = {
    "n_total": 239_695,
    "n_lum_limited": 62_220,
    "n_mass_limited": 41_801,
    "n_spirals_lum": 17_957,
    "arm_fractions_lum": {"1": 0.05, "2": 0.62, "3": 0.20, "4": 0.06, "5+": 0.06},
    "n_spirals_mass": 9_413,
    "arm_fractions_mass": {"1": 0.05, "2": 0.64, "3": 0.18, "4": 0.06, "5+": 0.07},
    "m2_fraction_lum": (0.621, 0.004),
    "m2_fraction_mass": (0.642, 0.004),
    "m2_fraction_highest_mass": (0.53, 0.02),     # log M* = 11.2 +/- 0.1
    "m5_fraction_mass": (0.068, 0.002),
    "m5_fraction_highest_mass": (0.15, 0.02),
    "m2_fraction_bluest": (0.32, 0.02),           # g - i = 0.67 +/- 0.07
    "m2_fraction_reddest": (0.84, 0.02),          # g - i = 1.17 +/- 0.05
    "mean_gi_by_arms_mass": {"1": 0.94, "2": 0.94, "3": 0.83, "4": 0.82, "5+": 0.82},
    "m2_fraction_densest": (0.75, 0.02),          # vs 0.643 +/- 0.005 overall
}

# Masters et al. 2011, MNRAS 411, 2026 (arXiv:1003.0449): bars in GZ2 disks.
# GZ2 data to July 2009; 0.01 < z < 0.06, M_r < -19.38 (H0 = 70, k-corrected to
# z = 0); log(a/b) < 0.3; >= 10 votes on the bar question; barred when
# N_bar >= N_nobar (raw votes).
MASTERS11 = {
    "n_disks": 13_665,
    "n_barred": 4_020,
    "bar_fraction": (0.294, 0.005),
    "bar_fraction_blue_disks": 0.20,    # (g - r) < 0.6, "constant ~20%"
    "bar_fraction_red_bulgy": 0.5,      # "over half of the red, bulge-dominated disks"
}

# Masters et al. 2010, MNRAS 405, 783 (arXiv:0910.4113): passive red spirals.
# GZ1, 0.03 < z < 0.085, M_r < -20.17; p_spiral >= 0.8 (debiased); face-on
# log(a/b) < 0.2; disky fracDeV <= 0.5; arms visible (GZ1 p_CW or p_ACW > 0.8);
# red: (g - r) > 0.63 - 0.02 (M_r + 20).
MASTERS10 = {
    "n_spirals": 5_433,
    "n_red": 294,
    "red_fraction": 0.06,
    "bar_fraction_red": (0.67, 0.05),   # "at least 67 +/- 5%" (70 +/- 5% in the abstract)
    "bar_fraction_blue": (0.27, 0.05),
    "red_line": (0.63, -0.02),          # (g - r) = 0.63 - 0.02 (M_r + 20)
}

# Schawinski et al. 2009, MNRAS 396, 818 (arXiv:0903.3415): blue early types.
# GZ1 clean early types (> 80% of weighted votes), 0.02 < z < 0.05,
# M_r < -20.7; blue = more than 3 sigma below a Gaussian fit to the red
# sequence in u - r, in 0.33-mag bins.
SCHAWINSKI09 = {
    "n_blue": 204,
    "blue_fraction": (0.057, 0.004),
}

# Bamford et al. 2009, MNRAS 393, 1324 (arXiv:0805.2612): morphology, colour,
# environment and stellar mass. Qualitative headline used for grading.
BAMFORD09 = {
    "headline": ("at fixed stellar mass, colour is highly sensitive to environment, "
                 "while morphology displays much weaker environmental trends"),
    "colour_density_driver": ("the colour-density relation is primarily driven by colour "
                              "at fixed morphology, especially red spirals at low mass"),
}

# Project D2 (07-cosmic-web), same SDSS DR17 Main Galaxy Sample.
D2 = {
    "color_divider": (-0.04996, -0.38606),   # rest-frame g - r = slope * M_r + intercept
    "red_fraction_sparsest_densest": (0.51, 0.81),
    "n_tracers": 134_975,
}
