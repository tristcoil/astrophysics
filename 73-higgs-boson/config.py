"""Project configuration: paths, the CMS open-data files, analysis settings, and
the measured Standard Model inputs our calculations start from.

Published results live here only for comparison at the end, never as inputs to
a fit. The exceptions are measured *inputs* -- particle masses, the Fermi
constant, the strong coupling -- which every calculation of the running
couplings needs, and the event weights of the CMS outreach analysis, which
define how its simulation is normalised (both labelled where they are used).
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"                 # the CMS ROOT files, byte for byte
CACHE_DIR = DATA_DIR / "cache"             # skims and histograms derived from them
LITERATURE_DIR = DATA_DIR / "literature"   # papers' source files and tables
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, CACHE_DIR, LITERATURE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-atlas-P7/1.0 (python urllib; educational use of CERN open data)"

# ---------------------------------------------------------------------------
# CMS open data (CERN Open Data portal, CC0)
# ---------------------------------------------------------------------------
#
# 2012 proton-proton collisions at sqrt(s) = 8 TeV, Run B and Run C, converted
# by Stefan Wunsch to a reduced NanoAOD format for education (records
# 12341 and 12361-12368, DOIs in the README). The portal serves the same EOS
# files over HTTPS with byte ranges; sizes and adler32 checksums are the ones
# the portal's record metadata lists, and every download is verified against
# them.

EOS_HTTPS = "https://opendata.cern.ch/eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/"

#: name -> (record id, sub-path on EOS, bytes, adler32, events)
CMS_FILES = {
    # Simulation (Monte Carlo) for the four-lepton analysis
    "SMHiggsToZZTo4L": (12361, "ForHiggsTo4Leptons/SMHiggsToZZTo4L.root",
                        42400229, "1dafc1d7", 299973),
    "ZZTo4mu": (12362, "ForHiggsTo4Leptons/ZZTo4mu.root", 184624480, "cb91f8d0", 1499064),
    "ZZTo4e": (12363, "ForHiggsTo4Leptons/ZZTo4e.root", 161660567, "f6fe4aff", 1499093),
    "ZZTo2e2mu": (12364, "ForHiggsTo4Leptons/ZZTo2e2mu.root", 173163793, "d03acd53", 1497445),
    # Collision data for the four-lepton analysis
    "Run2012B_DoubleMuParked": (12365, "ForHiggsTo4Leptons/Run2012B_DoubleMuParked.root",
                                3179316378, "daa50fff", 26084708),
    "Run2012C_DoubleMuParked": (12366, "ForHiggsTo4Leptons/Run2012C_DoubleMuParked.root",
                                4463307037, "2e91a64a", 35455705),
    "Run2012B_DoubleElectron": (12367, "ForHiggsTo4Leptons/Run2012B_DoubleElectron.root",
                                1817796893, "ab2abff8", 21474287),
    "Run2012C_DoubleElectron": (12368, "ForHiggsTo4Leptons/Run2012C_DoubleElectron.root",
                                2767819825, "31b461b8", 32537541),
    # Collision data for the dimuon spectrum (muon branches only, Run B + C)
    "Run2012BC_DoubleMuParked_Muons": (12341, "Run2012BC_DoubleMuParked_Muons.root",
                                       2244449133, "1fa61aca", 61540413),
}

SIMULATION = ("SMHiggsToZZTo4L", "ZZTo4mu", "ZZTo4e", "ZZTo2e2mu")
DATA_4L = ("Run2012B_DoubleMuParked", "Run2012C_DoubleMuParked",
           "Run2012B_DoubleElectron", "Run2012C_DoubleElectron")
DIMUON_FILE = "Run2012BC_DoubleMuParked_Muons"

#: The outreach analysis code (record 12360, Wunsch 2021, GPL-3.0), kept for reference.
OUTREACH_CODE_URL = ("https://opendata.cern.ch/eos/opendata/cms/software/"
                     "HiggsToFourLeptonsNanoAODOutreachAnalysis/")
OUTREACH_CODE_FILES = {"skim.cxx": "fcf2abdf", "histograms.py": "8caed224",
                       "plot.py": "ab61709f"}
DIMUON_CODE_URL = ("https://opendata.cern.ch/eos/opendata/cms/software/"
                   "DimuonSpectrumNanoAODOutreachAnalysis/DimuonSpectrumNanoAODOutreachAnalysis-1.4.zip")
DIMUON_CODE_ADLER32 = "5d4fc1b2"
#: The full package of record 12360, which includes the published plots as
#: vector PDFs: our answer key, bin by bin.
OUTREACH_ZIP_URL = OUTREACH_CODE_URL + "HiggsToFourLeptonsNanoAODOutreachAnalysis-1.1.zip"
OUTREACH_ZIP_ADLER32 = "334eb647"
OUTREACH_PLOTS_DIR = LITERATURE_DIR / "cms_outreach_code" / "plots"

#: Small literature files: name -> (URL, what for)
LITERATURE = {
    # Buttazzo et al. 2013, JHEP 12, 089: RGEs to three loops, matching at M_t,
    # effective potential, vacuum decay (formulas transcribed in rge.py)
    "buttazzo2013_src.tar.gz": ("https://arxiv.org/e-print/1307.3536v4", "RGEs, matching"),
    # Andreassen, Frost & Schwartz 2018, PRD 97, 056006 (arXiv v4): NLO vacuum decay rate
    "andreassen2018_src.tar.gz": ("https://arxiv.org/e-print/1707.08124v4", "decay-rate benchmark"),
    # D'Onofrio & Rummukainen 2016, PRD 93, 025003: the Higgs condensate across
    # the electroweak crossover from lattice simulations (the authors' data file)
    "donofrio_rummukainen2016_sm_data.txt": ("https://www.mv.helsinki.fi/home/rummukai/sm_data.txt",
                                            "Higgs condensate vs temperature"),
    # ATLAS Collaboration 2022, Nature 607, 52 (HEPData ins2104706, Figure 5):
    # measured Higgs couplings to W, Z, t, b, tau, mu
    "atlas2022_reduced_couplings.csv": (
        "https://www.hepdata.net/download/table/ins2104706/"
        "Figure%205%20without%20%24%5Ckappa_%7Bc%7D%24/1/csv", "couplings vs mass"),
    # SMDR 1.3 (Martin & Robertson 2019, arXiv:1907.02500): independent matching and running
    "smdr-1.3.tar.gz": ("https://www.niu.edu/spmartin/smdr/smdr-1.3.tar.gz", "cross-check code"),
}
ARXIV_USER_AGENT = "astrophysics-atlas/1.0 (educational; python urllib)"

#: Parallel byte-range download: chunk size and number of connections.
CHUNK_BYTES = 64 * 1024 * 1024
DOWNLOAD_WORKERS = 12

# ---------------------------------------------------------------------------
# The four-lepton analysis as published (record 12360, skim.cxx)
# ---------------------------------------------------------------------------
#
# Event weights put the simulation on the data's scale: integrated luminosity
# x cross-section / number of simulated events, with a 1.386 "correction of
# the simulation" for the ZZ background. These are the outreach analysis's own
# numbers, used as published.

LUMINOSITY_PB = 11.58 * 1000.0            # 11.58 fb^-1 of Run B + C
ZZ_SCALE_FACTOR = 1.386
XSEC_PB = {"SMHiggsToZZTo4L": 0.0065, "ZZTo4mu": 0.077, "ZZTo4e": 0.077, "ZZTo2e2mu": 0.18}
N_GENERATED = {"SMHiggsToZZTo4L": 299973.0, "ZZTo4mu": 1499064.0, "ZZTo4e": 1499093.0,
               "ZZTo2e2mu": 1497445.0}

#: Which final states each file contributes (as in skim.cxx).
FINAL_STATES = {
    "SMHiggsToZZTo4L": ("4mu", "4e", "2e2mu"),
    "ZZTo4mu": ("4mu",), "ZZTo4e": ("4e",), "ZZTo2e2mu": ("2e2mu",),
    "Run2012B_DoubleMuParked": ("4mu", "2e2mu"), "Run2012C_DoubleMuParked": ("4mu", "2e2mu"),
    "Run2012B_DoubleElectron": ("4e", "2e2mu"), "Run2012C_DoubleElectron": ("4e", "2e2mu"),
}

Z_MASS_NOMINAL = 91.2                       # GeV, the pairing target in skim.cxx
M4L_HIST = (36, 70.0, 180.0)                # the published histogram: 36 bins, 70-180 GeV

#: Mass of the Higgs boson in the signal simulation (GluGluToHToZZTo4L_M-125).
MC_HIGGS_MASS = 125.0

#: Branches the lepton analyses read.
MUON_BRANCHES = ("nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
                 "Muon_pfRelIso04_all", "Muon_dxy", "Muon_dz", "Muon_dxyErr", "Muon_dzErr")
ELECTRON_BRANCHES = ("nElectron", "Electron_pt", "Electron_eta", "Electron_phi",
                     "Electron_mass", "Electron_charge", "Electron_pfRelIso03_all",
                     "Electron_dxy", "Electron_dz", "Electron_dxyErr", "Electron_dzErr")
EVENT_ID_BRANCHES = ("run", "luminosityBlock", "event")

# ---------------------------------------------------------------------------
# Measured Standard Model inputs: PDG 2026 (Review of Particle Physics 2026; the
# values of its Monte Carlo mass table, shipped with the `particle` package,
# and of its physical-constants table). The spec quoted PDG 2024
# (125.20 +- 0.11, 172.57 +- 0.29); both are run in stability.py.
# ---------------------------------------------------------------------------

M_HIGGS = (125.13, 0.11)                    # GeV
M_TOP = (172.60, 0.27)                      # GeV, direct measurements
M_W = (80.362, 0.008)                       # GeV
M_Z = (91.1879, 0.0020)                     # GeV
ALPHA_S_MZ = (0.1180, 0.0009)
G_FERMI = 1.1663785e-5                      # GeV^-2 (muon lifetime)
M_PLANCK = 1.220890e19                      # GeV, sqrt(hbar c / G)
REDUCED_PLANCK = 2.435323e18                # GeV, M_Pl / sqrt(8 pi)
HBAR_GEV_S = 6.582119569e-25                # GeV s
