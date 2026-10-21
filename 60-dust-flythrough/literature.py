"""Published values this project checks against, each with its source."""

#: Edenhofer et al. (2024, A&A 685, A82): their Fig. 4 colour bar, read from the arXiv source image
#: (validation.fig4_check fits the tick marks), and statements in the text.
EDENHOFER = {
    "fig4_faceon_vmax_av": 1.698,          # 99.9% quantile of A_V through |z| <= 400 pc (our reading)
    "radial_resolution_pc": (0.4, 7.0),    # "highest/lowest distance discretization"
    "std_over_mean_projection": 0.10,      # "on the order of 10% of the posterior mean"
    "a_v_per_zgr23": 2.8,
}

#: O'Neill et al. (2024, ApJ 973, 136), Table 1 and section 3.
ONEILL = {
    "d_peak_pc": {"min": 71, "p2.28": 91, "median": 170, "p97.72": 362, "max": 616},
    "thickness_pc": {"min": 3, "p2.28": 13, "median": 35, "p97.72": 123, "max": 357},
    "n_peak_cm3": {"min": 0.02, "p2.28": 0.04, "median": 0.61, "p97.72": 21, "max": 770},
    "mass_msun": 6.0e5, "mass_sys_msun": 0.7e5, "mass_a09_msun": 8.0e5,
    "interior_volume_pc3": 1.9e7, "interior_equivalent_radius_pc": 165,
    "shell_extent_pc": {"x": (-300, 330), "y": (-355, 445), "z": (-300, 600)},
}

#: Zucker et al. (2022, Nature 601, 334), Methods (medians, 16th-84th percentiles).
ZUCKER2022 = {
    "t_exp_myr": (-14.39, -0.74, +0.78), "n0_cm3": (2.71, -1.02, +1.57), "dt_sne_myr": (1.06, -0.39, +0.63),
    "dR_pc": (23.31, -2.29, +2.54), "R_now_pc": (165, -6, +6), "v_now_kms": (6.7, -0.4, +0.5),
    "centre_pc": (39, 7, -18), "n_sne": (15, -7, +11), "shell_mass_msun": 1.4e6,
}

#: Alves et al. (2020, Nature 578, 237), Extended Data Tables 2 and 3.
ALVES2020 = {
    "length_kpc": (2.7, 0.2), "amplitude_pc": (160, 30), "scatter_pc": (60, 15), "mass_msun": 3e6,
    "average_period_kpc": 2.0,
    "params": {"x0": -910, "y0": -860, "z0": -30, "x1": -270, "y1": -20, "z1": -10, "x2": 290, "y2": 1400,
               "z2": 30, "P": 3560, "A": 160, "phi": 2.89, "gamma": 0.50, "delta": 6.68, "sigma": 62, "f": 0.22},
}

#: Konietzka et al. (2024, Nature 628, 62): text and Extended Data Table 2 (clouds-only fit).
KONIETZKA2024 = {
    "max_amplitude_pc": 220, "mean_wavelength_pc": 2000, "wavelength_range_pc": (400, 4000),
    "radius_pc": 47, "max_vertical_velocity_kms": 14, "radial_drift_kms": 5,
}

#: Zucker et al. (2020, A&A 633, A51): distances carry a ~5% systematic uncertainty.
ZUCKER2020 = {"systematic_fraction": 0.05, "n_sightlines": 326}

#: Hunt & Reffert (2023, A&A 673, A114).
HUNT_REFFERT = {"n_clusters": 7167, "n_reliable_arxiv": 4105, "n_reliable_abstract": 4114}
