"""Published nickel masses and related numbers we grade our results against.

These are never inputs to a fit. Values were checked against the papers'
own text or tables (arXiv sources), not secondary summaries. Distances differ
between papers; where it matters, the distance each paper assumed is noted
so the comparison can be made on equal terms.
"""

from __future__ import annotations

NICKEL = {
    "SN2011fe": [
        {"ref": "Pereira et al. 2013 (A&A 554, A27)", "m_ni": 0.53, "err": 0.11,
         "method": "Arnett's rule, alpha = 1.2", "dist_mpc": 6.4,
         "note": "peak L = (1.17 +/- 0.04) x 10^43 erg/s at 6.4 Mpc"},
    ],
    "SN2014J": [
        {"ref": "Churazov et al. 2014 (Nature 512, 406)", "m_ni": 0.61, "err": 0.13,
         "method": "INTEGRAL 56Co lines, days 50-100", "dist_mpc": 3.53},
        {"ref": "Diehl et al. 2015 (A&A 574, A72)", "m_ni": 0.49, "err": 0.09,
         "method": "INTEGRAL 56Co line light curves vs models", "dist_mpc": 3.3},
    ],
    "SN1987A": [
        {"ref": "Seitenzahl et al. 2014 (ApJ 792, 10)", "m_ni": 0.071, "err": 0.003,
         "method": "UVOIR light curve to 4334 days", "dist_mpc": 0.050},
    ],
    "SN2023ixf": [
        {"ref": "Hsu et al. 2025 (ApJ 990, 148)", "m_ni": 0.068, "err": 0.001,
         "method": "tail fit days 90-150, T0 = 190 +/- 3 d", "dist_mpc": 6.71},
        {"ref": "Singh et al. 2024 (ApJ 975, 132)", "m_ni": 0.059, "err": 0.001,
         "method": "tail fit days 100-150, T0 = 220 +/- 3 d", "dist_mpc": 6.82},
        {"ref": "Zimmerman et al. 2024 (Nature 627, 759)", "m_ni": 0.071, "err": None,
         "method": "as quoted by Hsu et al. 2025", "dist_mpc": None},
        {"ref": "Bersten et al. 2024 (A&A 681, L18)", "m_ni": 0.05, "err": None,
         "method": "hydrodynamical model, as quoted by Hsu et al. 2025", "dist_mpc": None},
        {"ref": "Moriya & Singh 2024 (PASJ 76, 1050)", "m_ni": 0.04, "err": None,
         "method": "model grid, as quoted by Hsu et al. 2025", "dist_mpc": None},
    ],
}

#: Gamma-ray escape times T0 measured from tails (days).
T0 = {
    "SN2023ixf": [("Hsu et al. 2025", 190.0, 3.0), ("Singh et al. 2024", 220.0, 3.0)],
}

#: The spec's headline checks.
CHECKS = {
    "SN2014J_gamma": (0.6, 0.1),  # Churazov et al. 2014, "0.6 +/- 0.1"
    "SN1987A": (0.07, None),
}
