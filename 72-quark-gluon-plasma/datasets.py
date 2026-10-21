"""The measured yields, compiled from HEPData and the papers themselves.

Each yield is dN/dy at mid-rapidity (|y| < 0.5): the number of particles per
collision in one unit of rapidity, with its statistical and systematic
uncertainty. For the fit the two are added in quadrature, as Andronic et al.
2018 did.

ALICE, Pb-Pb at 2.76 TeV, 0-10% most central collisions:

* pi, K, p, K0S, Lambda and phi were measured in 0-5% and 5-10%: the 0-10%
  value is their mean, with statistical errors combined in quadrature and
  systematic errors averaged (treated as fully correlated between the two).
* Xi, Omega, the deuteron and the hypertriton were measured in 0-10% itself.
* 3He was measured in 0-20%; like ALICE (and Thermal-FIST's compilation) we
  scale it to 0-10% with the deuteron's 0-10% / 0-20% ratio. The
  antideuteron and anti-3He come from the measured antinucleus/nucleus ratios.
* The hypertriton is published as yield x branching ratio into 3He pi-;
  ALICE's assumed BR of 25% converts it.
* Lambda-bar was not published separately: ALICE states that Lambda-bar/K0S
  equals Lambda/K0S. It is therefore not an independent measurement and is
  left out of our default fit (Andronic et al. plot it; a fit variant adds it).
* K*(892)0 is shown but never fitted: it lives only 4 fm/c, and pions
  scattering on its decay products after chemical freeze-out hide some of it.

STAR, Au+Au at 7.7-39 GeV, 0-5%: pi, K, p (whose yields are inclusive, i.e.
include protons from Lambda decays) and Lambda, Xi (corrected for weak-decay
feed-down), from the two beam-energy-scan papers.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from functools import lru_cache

import yaml

import config as cfg
import fetch_data

# Particle codes
PI, K, P, K0, LAM, XI, OM, PHI, KSTAR0 = 211, 321, 2212, 311, 3122, 3312, 3334, 333, 313
D, HE3, H3L, HE4, TRITON = 1000010020, 1000020030, 1010010030, 1000020040, 1000010030


@dataclass
class Yield:
    key: str
    label: str                 # matplotlib mathtext
    html: str                  # HTML label
    targets: dict              # {pdg: coefficient}
    feeddown: str              # "strong" or "weak"
    value: float
    stat: float
    sys: float
    fitted: bool = True
    group: str = ""            # "meson", "baryon", "nucleus", "resonance"
    source: str = ""
    note: str = ""

    @property
    def err(self) -> float:
        return math.hypot(self.stat, self.sys)

    @property
    def mass_pdg(self) -> int:
        return abs(next(iter(self.targets)))

    @property
    def is_anti(self) -> bool:
        return next(iter(self.targets)) < 0 and self.key not in ("pi-", "K-")


@dataclass
class Dataset:
    name: str
    title: str
    system: str               # "Pb" or "Au"
    sqrt_s: float             # GeV per nucleon pair
    centrality: str
    yields: list = field(default_factory=list)
    mu_q: str = "Q/B"         # "Q/B": from the nuclei's Q/B; "zero": mu_Q = 0
    notes: list = field(default_factory=list)
    charm_pairs: float = 0.0  # c-cbar pairs per unit y (their decays add pi and K)

    @property
    def q_over_b(self) -> float | None:
        return cfg.Q_OVER_B[self.system] if self.mu_q == "Q/B" else None

    def get(self, key: str) -> Yield:
        for y in self.yields:
            if y.key == key:
                return y
        raise KeyError(key)

    def fitted(self) -> list:
        return [y for y in self.yields if y.fitted]


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def hepdata_table(record: str, filename: str) -> dict:
    return yaml.safe_load((fetch_data.hepdata_dir(record) / filename).read_text())


def _errors(val: float, errs: list | None) -> tuple:
    """(stat, sys) from a HEPData error list; percentages resolved, several
    systematic components added in quadrature; an unlabelled error is total."""
    stat, sys2 = 0.0, 0.0
    for e in errs or []:
        if "symerror" in e:
            s = e["symerror"]
        else:
            s = 0.5 * (abs(e["asymerror"]["plus"]) + abs(e["asymerror"]["minus"]))
        if isinstance(s, str):
            s = float(s.strip().rstrip("%")) / 100.0 * val
        label = e.get("label", "")
        if label.startswith("stat"):
            stat = math.hypot(stat, s)
        else:
            sys2 += s * s
    return stat, math.sqrt(sys2)


def centrality_rows(record: str, filename: str, column: int = 0, factor: float = 1.0) -> dict:
    """{(lo, hi): (value, stat, sys)} for a table whose x is a centrality class."""
    t = hepdata_table(record, filename)
    iv = t["independent_variables"][0]["values"]
    dv = t["dependent_variables"][column]["values"]
    out = {}
    for x, v in zip(iv, dv):
        val = float(v["value"])
        st, sy = _errors(val, v.get("errors"))
        out[(float(x["low"]), float(x["high"]))] = (val * factor, st * factor, sy * factor)
    return out


def _combine(rows: dict, classes: list) -> tuple:
    """Mean of equal-width centrality classes: stat in quadrature, sys averaged."""
    vals = [rows[c] for c in classes]
    n = len(vals)
    v = sum(x[0] for x in vals) / n
    st = math.sqrt(sum(x[1] ** 2 for x in vals)) / n
    sy = sum(x[2] for x in vals) / n
    return v, st, sy


def _cent(c: str) -> tuple:
    lo, hi = c.split("-")
    return float(lo), float(hi)


def _split_classes(rows: dict, cent: str) -> list:
    """The published classes that tile one requested class (equal widths)."""
    lo, hi = _cent(cent)
    if (lo, hi) in rows:
        return [(lo, hi)]
    inside = sorted(c for c in rows if c[0] >= lo and c[1] <= hi)
    width = {c[1] - c[0] for c in inside}
    if not inside or inside[0][0] != lo or inside[-1][1] != hi or len(width) != 1:
        raise KeyError(f"cannot build {cent} from {sorted(rows)}")
    return inside


def from_classes(rows: dict, cent: str) -> tuple:
    return _combine(rows, _split_classes(rows, cent))


_PM = r"\s*\\pm\s*"
_NUM = r"([0-9.]+)"


def alice_nuclei_paper() -> dict:
    """Deuteron dN/dy per centrality, 3He in 0-20% and 20-80%, and the
    antinucleus/nucleus ratios, from the TeX of Adam et al. 2016 (PRC 93, 024917)."""
    tex = (fetch_data.arxiv_dir("1506.08951") / "antinuclei.tex").read_text(errors="replace")
    out = {"d": {}, "He3": {}, "dbar/d": {}, "He3bar/He3": {}}
    pat = re.compile(r"^(d|\$\^\{3\}\$He)\s*\((\d+)-(\d+)\\%\)\s*&.*?\$\s*\(" + _NUM + _PM + _NUM
                     + _PM + _NUM + r"\)\s*\\times\s*10\^\{(-?\d+)\}", re.M)
    for m in pat.finditer(tex):
        who = "d" if m.group(1) == "d" else "He3"
        scale = 10.0 ** int(m.group(7))
        out[who][(float(m.group(2)), float(m.group(3)))] = tuple(
            float(m.group(i)) * scale for i in (4, 5, 6))
    block = tex[tex.index("Anti-nuclei/nuclei"):tex.index("tab.:RatioTable")]
    current = None
    for line in block.splitlines():
        if "bar{\\rm d}" in line:
            current = "dbar/d"
        elif "overline{\\rm He}" in line:
            current = "He3bar/He3"
        m = re.search(r"(\d+)-(\d+)\\%\s*&\s*" + _NUM + r"\s*\$\\pm\$\s*" + _NUM
                      + r"\s*\$\\pm\$\s*" + _NUM, line)
        if m and current:
            out[current][(float(m.group(1)), float(m.group(2)))] = tuple(
                float(m.group(i)) for i in (3, 4, 5))
    if len(out["d"]) != 5 or len(out["He3"]) != 2 or len(out["dbar/d"]) != 5:
        raise RuntimeError(f"nuclei table parse failed: {out}")
    return out


def alice_he4_paper() -> dict:
    """(value, stat, sys) for 4He and anti-4He in 0-10%, from the abstract of
    Acharya et al. 2018 (NPA 971, 1); measured in |y| < 1, per unit rapidity."""
    tex = (fetch_data.arxiv_dir("1710.07531") / "alpha_paper.tex").read_text(errors="replace")
    pat = re.compile(r"\(" + _NUM + _PM + _NUM + r"\s*~?\(\\mathrm\{stat\}\)" + _PM + _NUM
                     + r"~?\(\\mathrm\{syst\}\)\)\\times\s*10\^\{(-\d+)\}")
    found = pat.findall(tex)
    if len(found) < 2:
        raise RuntimeError("4He yields not found")
    vals = [tuple(float(x) * 10.0 ** int(f[3]) for x in f[:3]) for f in found[:2]]
    return {"He4": vals[0], "He4bar": vals[1]}


# ---------------------------------------------------------------------------
# ALICE, Pb-Pb 2.76 TeV
# ---------------------------------------------------------------------------

ALICE_LABELS = {
    "pi+": (r"$\pi^+$", "π⁺"), "pi-": (r"$\pi^-$", "π⁻"),
    "K+": (r"K$^+$", "K⁺"), "K-": (r"K$^-$", "K⁻"), "K0S": (r"K$^0_S$", "K⁰<sub>S</sub>"),
    "phi": (r"$\phi$", "φ"), "p": (r"p", "p"), "pbar": (r"$\overline{\rm p}$", "p̄"),
    "Lambda": (r"$\Lambda$", "Λ"), "Lambdabar": (r"$\overline{\Lambda}$", "Λ̄"),
    "Xi-": (r"$\Xi^-$", "Ξ⁻"), "Xibar+": (r"$\overline{\Xi}^{+}$", "Ξ̄⁺"),
    "Omega-": (r"$\Omega^-$", "Ω⁻"), "Omegabar+": (r"$\overline{\Omega}^{+}$", "Ω̄⁺"),
    "d": (r"d", "d"), "dbar": (r"$\overline{\rm d}$", "d̄"),
    "He3": (r"$^3$He", "³He"), "He3bar": (r"$^3\overline{\rm He}$", "³He̅"),
    "H3L": (r"$^3_\Lambda$H", "³<sub>Λ</sub>H"),
    "H3Lbar": (r"$^3_{\overline{\Lambda}}\overline{\rm H}$", "³<sub>Λ̄</sub>H̅"),
    "He4": (r"$^4$He", "⁴He"), "He4bar": (r"$^4\overline{\rm He}$", "⁴He̅"),
    "K*0": (r"K$^{*0}$", "K*⁰"),
}


def _y(key, targets, conv, vss, group, source, fitted=True, note=""):
    lab, html = ALICE_LABELS[key]
    return Yield(key, lab, html, targets, conv, *vss, fitted=fitted, group=group,
                 source=source, note=note)


@lru_cache(maxsize=None)
def alice(cent: str = "0-10", include_lambdabar: bool = False) -> Dataset:
    """ALICE Pb-Pb 2.76 TeV yields in one centrality class."""
    lo, hi = _cent(cent)
    ys = []
    pikp = {name: centrality_rows("alice_pikp", "Table_31.yaml", col)
            for col, name in enumerate(["dNch", "pi+", "pi-", "K+", "K-", "p", "pbar"])}
    src = "HEPData ins1222333 Table 31"
    for key, pdg, grp in (("pi+", PI, "meson"), ("pi-", -PI, "meson"), ("K+", K, "meson"),
                          ("K-", -K, "meson"), ("p", P, "baryon"), ("pbar", -P, "baryon")):
        ys.append(_y(key, {pdg: 1.0}, "strong", from_classes(pikp[key], cent), grp, src))
    k0s = centrality_rows("alice_k0s_lambda", "Table22.yaml")
    ys.append(_y("K0S", {K0: 0.5, -K0: 0.5}, "strong", from_classes(k0s, cent), "meson",
                 "HEPData ins1243863 Table 22"))
    phi = centrality_rows("alice_kstar_phi", "Table21.yaml")
    ys.append(_y("phi", {PHI: 1.0}, "strong", from_classes(phi, cent), "meson",
                 "HEPData ins1288320 Table 21"))
    lam = centrality_rows("alice_k0s_lambda", "Table23.yaml")
    lam_v = from_classes(lam, cent)
    ys.append(_y("Lambda", {LAM: 1.0}, "strong", lam_v, "baryon", "HEPData ins1243863 Table 23",
                 note="includes Sigma0 -> Lambda gamma; Xi feed-down subtracted"))
    ys.append(_y("Lambdabar", {-LAM: 1.0}, "strong", lam_v, "baryon",
                 "= Lambda (ALICE: Lambda-bar/K0S = Lambda/K0S)", fitted=include_lambdabar,
                 note="not published separately; not an independent measurement"))
    xo = {name: centrality_rows("alice_xi_omega", "Table11.yaml", col)
          for col, name in ((0, "Xi-"), (1, "Xibar+"), (3, "Omega-"), (4, "Omegabar+"))}
    for key, pdg in (("Xi-", XI), ("Xibar+", -XI), ("Omega-", OM), ("Omegabar+", -OM)):
        ys.append(_y(key, {pdg: 1.0}, "strong", from_classes(xo[key], cent), "baryon",
                     "HEPData ins1243865 Table 11"))
    nuc = alice_nuclei_paper()
    src_n = "Adam et al. 2016 (PRC 93, 024917) Tables III-IV"
    d = from_classes(nuc["d"], cent) if (lo, hi) not in nuc["d"] else nuc["d"][(lo, hi)]
    ys.append(_y("d", {D: 1.0}, "strong", d, "nucleus", src_n))
    r = nuc["dbar/d"][(lo, hi)] if (lo, hi) in nuc["dbar/d"] else None
    if r is not None:
        rel = math.hypot(d[2] / d[0], r[2] / r[0])
        ys.append(_y("dbar", {-D: 1.0}, "strong",
                     (d[0] * r[0], d[0] * r[0] * math.hypot(d[1] / d[0], r[1] / r[0]),
                      d[0] * r[0] * rel), "nucleus", src_n + ", d x (dbar/d)"))
    if cent == "0-10":
        # 3He: 0-20% scaled to 0-10% with the deuteron's centrality dependence
        h = nuc["He3"][(0.0, 20.0)]
        d20 = (nuc["d"][(0.0, 10.0)][0] + nuc["d"][(10.0, 20.0)][0]) / 2
        s = nuc["d"][(0.0, 10.0)][0] / d20
        ys.append(_y("He3", {HE3: 1.0}, "strong", (h[0] * s, h[1] * s, h[2] * s), "nucleus",
                     src_n + ", 0-20% x d(0-10%)/d(0-20%)",
                     note=f"scaled from 0-20% by {s:.4f}"))
        r3 = nuc["He3bar/He3"][(0.0, 20.0)]
        v = h[0] * s * r3[0]
        ys.append(_y("He3bar", {-HE3: 1.0}, "strong",
                     (v, v * math.hypot(h[1] / h[0], r3[1] / r3[0]),
                      v * math.hypot(h[2] / h[0], r3[2] / r3[0])), "nucleus",
                     src_n + ", 3He x (anti-3He/3He)"))
        br = cfg.HYPERTRITON_BR_3HE_PI
        for col, key, pdg in ((1, "H3L", H3L), (2, "H3Lbar", -H3L)):
            rows = centrality_rows("alice_hypertriton", "Table4.yaml", col, factor=1e-5 / br)
            ys.append(_y(key, {pdg: 1.0}, "strong", rows[(0.0, 10.0)], "nucleus",
                         "HEPData ins1380234 Table 4 / BR(25%)"))
        he4 = alice_he4_paper()
        src4 = "Acharya et al. 2018 (NPA 971, 1), abstract"
        ys.append(_y("He4", {HE4: 1.0}, "strong", he4["He4"], "nucleus", src4))
        ys.append(_y("He4bar", {-HE4: 1.0}, "strong", he4["He4bar"], "nucleus", src4))
        # K*0: measured in 0-20%; carried to 0-10% with the measured K*/K- ratio
        kk = centrality_rows("alice_kstar_phi", "Table22.yaml")[(0.0, 20.0)]
        km = from_classes(pikp["K-"], "0-10")
        v = kk[0] * km[0]
        ys.append(_y("K*0", {KSTAR0: 0.5, -KSTAR0: 0.5}, "strong",
                     (v, v * math.hypot(kk[1] / kk[0], 0.0),
                      v * math.hypot(kk[2] / kk[0], km[2] / km[0])), "resonance",
                     "HEPData ins1288320 Table 22 (K*/K-, 0-20%) x K-(0-10%)", fitted=False,
                     note="short-lived (4 fm/c): depleted by rescattering"))
    return Dataset(f"alice_{cent}", f"ALICE Pb–Pb √s_NN = 2.76 TeV, {cent}%", "Pb", 2760.0,
                   cent, ys, charm_pairs=alice_charm_pairs(cent))


def alice_npart() -> dict:
    """<N_part> per centrality class (HEPData ins1243865 Table 11)."""
    t = hepdata_table("alice_xi_omega", "Table11.yaml")
    c = t["independent_variables"][0]["values"]
    n = t["independent_variables"][1]["values"]
    return {(float(a["low"]), float(a["high"])): float(b["value"]) for a, b in zip(c, n)}


def alice_charm_pairs(cent: str) -> float:
    """c-cbar pairs per unit y: cfg.CHARM_PAIRS_PER_DY in 0-10%, scaled to other
    classes like binary collisions, N_coll ~ N_part^(4/3) (Glauber-like)."""
    npart = alice_npart()
    ref = npart[(0.0, 10.0)]
    return cfg.CHARM_PAIRS_PER_DY[0] * (npart[_cent(cent)] / ref) ** (4 / 3)


ALICE_CENTRALITIES = ("0-10", "10-20", "20-40", "40-60", "60-80")


@lru_cache(maxsize=None)
def alice_andronic_compilation() -> Dataset:
    """The 0-10% yields exactly as Andronic et al. 2018 fitted them.

    Their Fig. 1 (read back by figparse) shows the values of Thermal-FIST's
    ALICE 0-10% file to four digits, including the paper-rounded K0S (100)
    and p (31.5) and a separate Lambda-bar (24.27). Used only to separate
    data-compilation differences from model differences.
    """
    import copy

    tf = thermal_fist_compilation()
    base = alice("0-10", include_lambdabar=True)
    ys = []
    for y in base.yields:
        y2 = copy.copy(y)
        pdg = 310 if y.key == "K0S" else next(iter(y.targets))
        if pdg in tf:
            y2.value, y2.stat, y2.sys = tf[pdg][0], 0.0, tf[pdg][1]
            y2.source = "Thermal-FIST ALICE 0-10% file (= Andronic et al. 2018 Fig. 1)"
        ys.append(y2)
    return Dataset("alice_0-10_andronic", "ALICE 0-10% as compiled by Andronic et al. 2018",
                   "Pb", 2760.0, "0-10", ys, charm_pairs=base.charm_pairs)


def thermal_fist_compilation() -> dict:
    """Thermal-FIST's ALICE 0-10% file (cross-check only): {pdg: (value, error)}."""
    path = cfg.LISTS_DIR / "input/data/ALICE-PbPb2.76TeV-0-10-all.dat"
    out = {}
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("#") or not line.strip():
            continue
        p = line.split()
        out[int(p[1])] = (float(p[5]), float(p[6]))
    return out


# ---------------------------------------------------------------------------
# STAR, Au+Au 7.7-39 GeV (beam energy scan), 0-5%
# ---------------------------------------------------------------------------

STAR_ENERGIES = (7.7, 11.5, 19.6, 27.0, 39.0)
_STAR_FILE = {7.7: "7", 11.5: "11", 19.6: "19", 27.0: "27", 39.0: "39"}


def star_npart_0_5() -> dict:
    """<Npart> of the 0-5% class at each energy (STAR 2017, Table II)."""
    t = hepdata_table("star_bes_bulk", "table2.yaml")
    out = {}
    for dv in t["dependent_variables"]:
        quals = {q["name"]: q["value"] for q in dv.get("qualifiers", [])}
        energy = None
        for v in quals.values():
            m = re.search(r"([0-9.]+)\s*GeV", str(v))
            if m:
                energy = float(m.group(1))
        name = dv["header"]["name"]
        m = re.search(r"([0-9.]+)\s*GeV", name)
        if energy is None and m:
            energy = float(m.group(1))
        if energy is None:
            continue
        val = dv["values"][0]
        out[energy] = (float(val["value"]), _errors(float(val["value"]), val.get("errors"))[1])
    return out


@lru_cache(maxsize=None)
def star(energy: float) -> Dataset:
    """STAR 0-5% Au+Au: pi, K, p (|y| < 0.1), Lambda, Xi (|y| < 0.5)."""
    npart = star_npart_0_5()[energy][0]
    t = hepdata_table("star_bes_bulk", "fig21.yaml")
    ys = []
    names = {"pi-": (-PI, "strong"), "pi+": (PI, "strong"), "ka-": (-K, "strong"),
             "ka+": (K, "strong"), "pbar": (-P, "weak"), "proton": (P, "weak")}
    keymap = {"pi-": "pi-", "pi+": "pi+", "ka-": "K-", "ka+": "K+", "pbar": "pbar", "proton": "p"}
    for dv in t["dependent_variables"]:
        who = dv["qualifiers"][0]["value"]
        pdg, conv = names[who]
        for x, v in zip(t["independent_variables"][0]["values"], dv["values"]):
            if abs(float(x["value"]) - energy) < 1e-6:
                val = float(v["value"])
                st, sy = _errors(val, v.get("errors"))
                f = npart / 2
                key = keymap[who]
                lab, html = ALICE_LABELS[key]
                note = "inclusive: includes weak-decay feed-down" if conv == "weak" else ""
                ys.append(Yield(key, lab, html, {pdg: 1.0}, conv, val * f, st * f, sy * f,
                                group="baryon" if "p" in key and key != "pi+" and key != "pi-"
                                else "meson",
                                source="HEPData ins1510593 fig. 21 x <Npart>/2", note=note))
    files = {"Lambda": ("lambda", LAM), "Lambdabar": ("lambdabar", -LAM),
             "Xi-": ("xi", XI), "Xibar+": ("xibar", -XI)}
    for key, (stem, pdg) in files.items():
        t2 = hepdata_table("star_bes_strange", f"fig17_{stem}_auau{_STAR_FILE[energy]}.yaml")
        iv = [float(x["value"]) for x in t2["independent_variables"][0]["values"]]
        j = max(range(len(iv)), key=lambda i: iv[i])          # most central = 0-5%
        v = t2["dependent_variables"][0]["values"][j]
        val = float(v["value"])
        st, sy = _errors(val, v.get("errors"))
        lab, html = ALICE_LABELS[key]
        ys.append(Yield(key, lab, html, {pdg: 1.0}, "strong", val, st, sy, group="baryon",
                        source="HEPData ins1738953 fig. 17", note="weak-decay feed-down subtracted"))
    return Dataset(f"star_{energy:g}", f"STAR Au+Au √s_NN = {energy:g} GeV, 0-5%", "Au",
                   energy, "0-5", ys)


def star_gce_inputs(energy: float) -> tuple:
    """STAR's own GCE fit (2017, Fig. 25): the yields it used and its model values."""
    t = hepdata_table("star_bes_bulk", "fig25_gcey.yaml")
    parts = [x["value"] for x in t["independent_variables"][0]["values"]]
    data, model = {}, {}
    for dv in t["dependent_variables"]:
        quals = {q["name"]: q["value"] for q in dv["qualifiers"]}
        if quals.get("Energy") != f"{energy:g} GeV":
            continue
        kind = dv["header"]["name"]
        for p, v in zip(parts, dv["values"]):
            if "data" in kind:
                data[p] = (float(v["value"]), float(v["errors"][0]["symerror"]))
            elif "model" in kind:
                model[p] = float(v["value"])
    return data, model


@lru_cache(maxsize=None)
def star_gce_dataset(energy: float) -> Dataset:
    """The exact inputs of STAR's own grand-canonical fit (2017, Fig. 25),
    with STAR's conventions: inclusive protons, mu_Q = 0."""
    data, _ = star_gce_inputs(energy)
    keys = {"pi+": ("pi+", PI, "strong"), "pi-": ("pi-", -PI, "strong"),
            "ka+": ("K+", K, "strong"), "ka-": ("K-", -K, "strong"),
            "p": ("p", P, "weak"), "pbar": ("pbar", -P, "weak"),
            "lambda": ("Lambda", LAM, "strong"), "lambdabar": ("Lambdabar", -LAM, "strong"),
            "Xi": ("Xi-", XI, "strong"), "Xibar": ("Xibar+", -XI, "strong")}
    ys = []
    for name, (val, err) in data.items():
        key, pdg, conv = keys[name]
        lab, html = ALICE_LABELS[key]
        ys.append(Yield(key, lab, html, {pdg: 1.0}, conv, val, 0.0, err,
                        group="meson" if key[0] in "pK" and key not in ("p", "pbar") else "baryon",
                        source="HEPData ins1510593 fig. 25 (STAR's GCE fit inputs)"))
    return Dataset(f"star_gce_{energy:g}", f"STAR GCE inputs, {energy:g} GeV", "Au", energy,
                   "0-5", ys, mu_q="zero")


def star_published_freezeout() -> dict:
    """STAR's GCE (T_ch, mu_B) for 0-5% central collisions (2017, Fig. 33a)."""
    t = hepdata_table("star_bes_bulk", "fig33_0-5%.yaml")
    e = [float(x["value"]) for x in t["independent_variables"][0]["values"]]
    out = {}
    cols = {dv["header"]["name"]: dv["values"] for dv in t["dependent_variables"]}
    mub = cols["muB [MeV]"]
    tch = cols["Tch [MeV]"]
    for i, en in enumerate(e):
        out[en] = {"muB": float(mub[i]["value"]), "muB_err": float(mub[i]["errors"][0]["symerror"]),
                   "T": float(tch[i]["value"]), "T_err": float(tch[i]["errors"][0]["symerror"])}
    return out
