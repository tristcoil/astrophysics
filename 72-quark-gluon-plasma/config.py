"""Project configuration: paths, data sources, physical constants, model
settings, and the published values our results are graded against.

Published numbers live here only for comparison at the end, never as inputs to
a fit. The exceptions are labelled where they are used:

* measured hadron yields (they are the data),
* the hypertriton's 25% branching ratio into 3He + pi-, which ALICE assumed
  to turn its measured yield x BR into a yield,
* lattice QCD results (the equation of state, T_pc and the curvature of the
  crossover line): they are what the freeze-out temperature is compared with,
* the baryon density of the universe (Planck's omega_b), which sets the
  early universe's tiny baryon chemical potential.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
DERIVED_DIR = DATA_DIR / "derived"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

HEPDATA_DIR = DATA_DIR / "hepdata"
ARXIV_DIR = DATA_DIR / "arxiv"
LISTS_DIR = DATA_DIR / "thermal_fist"

for _d in (DATA_DIR, FIGURE_DIR, DERIVED_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: HEPData sits behind Cloudflare, which challenges browser user agents and
#: refuses Python's default one; this honest string is let through by all three
#: hosts we use (HEPData, arXiv, GitHub).
USER_AGENT = "python-requests/2.32.3 (astrophysics-atlas)"

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: HEPData records (INSPIRE id, record version), downloaded as the full YAML
#: submission (a few kB each). Versions are pinned so numbers cannot drift.
HEPDATA = {
    "alice_pikp": (1222333, 1,
                   "ALICE: pi, K, p in Pb-Pb at 2.76 TeV vs centrality "
                   "(Abelev et al. 2013, PRC 88, 044910)"),
    "alice_k0s_lambda": (1243863, 1,
                         "ALICE: K0S and Lambda in Pb-Pb at 2.76 TeV "
                         "(Abelev et al. 2013, PRL 111, 222301)"),
    "alice_xi_omega": (1243865, 1,
                       "ALICE: Xi and Omega in Pb-Pb at 2.76 TeV "
                       "(Abelev et al. 2014, PLB 728, 216)"),
    "alice_kstar_phi": (1288320, 1,
                        "ALICE: K*(892)0 and phi(1020) in Pb-Pb at 2.76 TeV "
                        "(Abelev et al. 2015, PRC 91, 024609)"),
    "alice_nuclei": (1380491, 2,
                     "ALICE: light nuclei in pp and Pb-Pb "
                     "(Adam et al. 2016, PRC 93, 024917)"),
    "alice_hypertriton": (1380234, 1,
                          "ALICE: hypertriton in Pb-Pb at 2.76 TeV "
                          "(Adam et al. 2016, PLB 754, 360)"),
    "alice_photons": (1394677, 1,
                      "ALICE: direct photons in Pb-Pb at 2.76 TeV "
                      "(Adam et al. 2016, PLB 754, 235)"),
    "star_bes_bulk": (1510593, 1,
                      "STAR: bulk properties from the Beam Energy Scan, "
                      "Au+Au 7.7-39 GeV (Adamczyk et al. 2017, PRC 96, 044904)"),
    "star_bes_strange": (1738953, 2,
                         "STAR: strange hadrons in Au+Au 7.7-39 GeV "
                         "(Adam et al. 2020, PRC 102, 034909)"),
}

#: arXiv e-print sources (TeX + vector figures). Used for tables that are not
#: on HEPData and for the published answer key.
ARXIV = {
    "1710.09425": "Andronic, Braun-Munzinger, Redlich & Stachel 2018, Nature 561, 321 "
                  "(the statistical hadronization fit: answer key, figures read back)",
    "1303.0737": "ALICE 2013, PRC 88, 044910 (pi, K, p yields table)",
    "1307.5530": "ALICE 2013, PRL 111, 222301 (K0S and Lambda table)",
    "1506.08951": "ALICE 2016, PRC 93, 024917 (deuteron and 3He yields, antinucleus ratios)",
    "1710.07531": "ALICE 2018, NPA 971, 1 (4He and anti-4He yields)",
    "1509.07324": "ALICE 2016, PLB 754, 235 (direct photons: published slopes)",
    "1309.5258": "Borsanyi et al. 2014, PLB 730, 99 (Wuppertal-Budapest equation of state)",
    "1407.6387": "HotQCD 2014, PRD 90, 094503 (equation of state, Table 1)",
    "1203.0784": "HotQCD 2012, PRD 86, 034509 (charge susceptibilities, Table V)",
    "1812.08235": "HotQCD 2019, PLB 795, 15 (chiral crossover: T_pc, curvature)",
    "2002.02821": "Borsanyi et al. 2020, PRL 125, 052001 (T_c(mu_B))",
    "1701.07065": "STAR 2017, PRC 96, 044904 (BES conventions: inclusive protons)",
}

#: Thermal-FIST (Vovchenko & Stoecker 2019, CPC 244, 295): we take only its
#: hadron lists and decay tables (text files), pinned to one commit. The ALICE
#: compilation and the lattice transcriptions are read only as cross-checks of
#: our own compilations from HEPData and the papers.
THERMAL_FIST_COMMIT = "4086191e952285cd572ba2bbbefb2fdc55401c9b"
THERMAL_FIST_RAW = ("https://raw.githubusercontent.com/vlvovch/Thermal-FIST/"
                    f"{THERMAL_FIST_COMMIT}/")
THERMAL_FIST_FILES = (
    "input/list/PDG2025/list-withnuclei.dat",
    "input/list/PDG2025/list-withcharm.dat",
    "input/list/PDG2025/decays.dat",
    "input/list/PDG2014/list-withnuclei.dat",
    "input/list/PDG2014/decays.dat",
    "input/data/ALICE-PbPb2.76TeV-0-10-all.dat",
    "input/lqcd/HotQCD-EoS.dat",
    "input/lqcd/WB-EoS.dat",
)
HADRON_LISTS = {
    "PDG2025": ("input/list/PDG2025/list-withnuclei.dat", "input/list/PDG2025/decays.dat"),
    "PDG2014": ("input/list/PDG2014/list-withnuclei.dat", "input/list/PDG2014/decays.dat"),
    # charm hadrons, only to count the pions and kaons their decays add
    "PDG2025charm": ("input/list/PDG2025/list-withcharm.dat", "input/list/PDG2025/decays.dat"),
}
#: Charm-anticharm pairs per unit rapidity in 0-10% Pb-Pb at 2.76 TeV, for the
#: charm feed-down variant: ALICE's total pp charm cross-section at 2.76 TeV,
#: 4.8 (+2.4/-1.7) mb (JHEP 07 (2012) 191), about 1/8 of it at mid-rapidity,
#: times T_AA = 23.4 /mb and ~0.75 for nuclear shadowing. A rough estimate,
#: hence the 50% error; Andronic et al. 2018 include this contribution.
CHARM_PAIRS_PER_DY = (10.0, 5.0)

# ---------------------------------------------------------------------------
# Physical constants (natural units: GeV, fm)
# ---------------------------------------------------------------------------

HBARC = 0.1973269804            # GeV fm
HBAR_GEV_S = 6.582119569e-25    # GeV s
M_PLANCK = 1.220890e19          # GeV (Planck mass, G = 1/M_Pl^2)
MEV_TO_K = 1.160451812e10       # 1 MeV / k_B in kelvin
#: Mass number and charge of the colliding nuclei: Q/B of the fireball.
Q_OVER_B = {"Pb": 82 / 208, "Au": 79 / 197}

# ---------------------------------------------------------------------------
# Model settings
# ---------------------------------------------------------------------------

DEFAULT_LIST = "PDG2025"
#: Resonance widths: relativistic Breit-Wigner over [max(threshold, m - 2G),
#: m + 2G], as Thermal-FIST's "BWTwoGamma" (its own examples' default); the
#: alternative is the zero-width approximation.
DEFAULT_WIDTHS = "bw"
#: Excluded-volume radius used by Andronic et al. 2018 for all hadrons; it
#: only rescales the volume (all densities fall by the same factor).
EXCLUDED_VOLUME_R0_FM = 0.3
#: Hypertriton -> 3He pi- branching ratio assumed by ALICE (Kamada et al. 1998)
#: to convert its measured yield x BR into a yield.
HYPERTRITON_BR_3HE_PI = 0.25

# ---------------------------------------------------------------------------
# The early universe
# ---------------------------------------------------------------------------

#: Planck 2018 (TT,TE,EE+lowE+lensing+BAO, Table 2): physical baryon density,
#: and the CMB temperature today.
PLANCK_OMEGA_B = (0.02242, 0.00014)
T_CMB_K = 2.7255

# ---------------------------------------------------------------------------
# Earlier atlas projects this one connects to
# ---------------------------------------------------------------------------

P2_DIR = ROOT_DIR / "68-thermal-history"   # P2: g*(T) and the cosmic calendar
P5_DIR = ROOT_DIR / "71-antimatter"        # P5: the baryon asymmetry
L4_DIR = ROOT_DIR / "29-element-origins"   # L4: where the elements come from

# ---------------------------------------------------------------------------
# Published values (answer key: compared with, never fitted)
# ---------------------------------------------------------------------------

#: Andronic, Braun-Munzinger, Redlich & Stachel 2018 (Nature 561, 321), text
#: and Fig. 3 caption, read in the arXiv source (1710.09425).
ANDRONIC_2018 = {
    "T": 156.5, "T_err": 1.5,             # MeV
    "muB": 0.7, "muB_err": 3.8,           # MeV
    "V": 5280.0, "V_err": 410.0,          # fm^3, with excluded volume r0 = 0.3 fm
    "proton_deviation_sigma": 2.7,
    "T_nuclei": 159.0, "T_nuclei_err": 5.0,
    "pion_decay_fraction": 0.70,          # share of pions from resonance decays
    "T_model_systematic_fraction": 0.03,  # "at most 3%"
    # Energy dependence: T = T_lim / (1 + exp(2.60 - ln(sqrt s)/0.45)),
    # mu_B = a / (1 + 0.288 sqrt s), sqrt s in GeV.
    "param_T_lim": 158.4, "param_T_lim_err": 1.4, "param_a": 1307.5,
}
#: HotQCD 2019 (Bazavov et al., PLB 795, 15), abstract and Table 2 (B,f column:
#: strangeness neutral, n_Q = 0.4 n_B, as at freeze-out).
HOTQCD_2019 = {"Tpc": 156.5, "Tpc_err": 1.5, "kappa2": 0.012, "kappa2_err": 0.004,
               "kappa4": 0.000, "kappa4_err": 0.004,
               "eps_pc": 0.42, "eps_pc_err": 0.06,   # GeV/fm^3 at T_pc
               "s_pc": 3.7, "s_pc_err": 0.5}         # fm^-3 at T_pc
#: Borsanyi et al. 2020 (PRL 125, 052001), strangeness neutral.
BORSANYI_2020 = {"Tc": 158.0, "Tc_err": 0.6, "kappa2": 0.0153, "kappa2_err": 0.0018,
                 "kappa4": 0.00032, "kappa4_err": 0.00067}
#: ALICE's own thermal-model statements.
ALICE_NUCLEI_2016 = {"T_p_d_3He": 156.0, "T_err": 4.0, "chi2_ndf": 0.4}
ALICE_HE4_2018 = {"T_nuclei": 154.0, "T_err": 4.0}
#: ALICE direct photons (PLB 754, 235), 0-20%, exponential fit 0.9 < pT < 2.1 GeV/c.
ALICE_PHOTONS_2016 = {"T_eff": 297.0, "stat": 12.0, "syst": 41.0,
                      "T_eff_no_pqcd_subtraction": 304.0, "stat_ns": 11.0, "syst_ns": 40.0,
                      "fit_range": (0.9, 2.1)}
