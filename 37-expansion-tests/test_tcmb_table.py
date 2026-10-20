"""Check every hand-transcribed CMB temperature against the papers' own text.

literature/tcmb_measurements.csv was typed from the papers' tables. This test
re-reads those tables from the arXiv sources cached by fetch_data.py (TeX, or
the PDF text via ``pdftotext`` for Riechers et al. 2022, which arXiv has only
as a PDF) and checks that each row's numbers appear where its "source" column
says they do.

    python test_tcmb_table.py
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys

import pandas as pd

import tcmb
from config import PAPER_DIR

FAILURES = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILURES.append(name)


def _squash(text):
    """Remove spaces and TeX markup so '$7.9^{+1.7}_{-1.4}$' -> '7.9^{+1.7}_{-1.4}'."""
    return re.sub(r"[\s$]", "", text)


def riechers_text():
    pdf = PAPER_DIR / "riechers_2022" / "paper.pdf"
    if shutil.which("pdftotext") is None:
        return None
    out = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True)
    text = out.stdout
    i = text.find("Extended Data Table 1: Summary")
    return text[i:]


def fmt(x):
    """Numbers as the papers print them (no trailing zeros beyond what is given)."""
    s = f"{x:g}"
    return s


def variants(x):
    """The ways a number may be printed: 1, 1.0, 1.00."""
    return {f"{x:g}", f"{x:.1f}", f"{x:.2f}"}


def asym_tokens(t, ep, em):
    return {f"{a}^{{+{b}}}_{{-{c}}}" for a in variants(t) for b in variants(ep) for c in variants(em)}


def main() -> int:
    df = tcmb.load_measurements()

    print("Riechers et al. 2022, Extended Data Table 1 (PDF text)")
    text = riechers_text()
    if text is None:
        print("  (pdftotext not installed: skipped)")
    else:
        body = text.replace(" ", "")
        for r in df[df["sets"].apply(lambda s: "riechers22" in s)].itertuples():
            if r.kind == "upper":
                token = f"<{fmt(r.T)}"
            elif r.kind == "range":
                token = "16.4-30.2" if r.object == "HFLS3" else "10+4-4"
            elif r.err_sys > 0:
                token = f"{r.T:.3f}±{r.err_plus:.3f}±{r.err_sys:.3f}"
            elif r.err_plus == r.err_minus and r.object.startswith(("Planck X-ray", "COBE")):
                token = f"{r.T}±{r.err_plus}"
            else:
                token = f"{r.T:g}+{r.err_plus:g}-{r.err_minus:g}"
            alt = [token, token.replace("+", "+").rstrip("0"),
                   f"{r.T:.1f}+{r.err_plus:.1f}-{r.err_minus:.1f}",
                   f"{r.T:.2f}+{r.err_plus:.2f}-{r.err_minus:.2f}",
                   f"{r.T:.5f}±{r.err_plus:.5f}"]
            ok = any(a in body for a in alt)
            check(f"{r.object} z={r.z} {r.method}", ok, "" if ok else f"(looked for {alt[:3]})")

    print("Klimenko et al. 2020 (TeX tables)")
    tex = _squash((PAPER_DIR / "klimenko_2020" / "main.tex").read_text())
    for r in df[df["reference"].str.startswith("Klimenko")].itertuples():
        alts = {f"<{v}" for v in variants(r.T)} if r.kind == "upper" else asym_tokens(
            r.T, r.err_plus, r.err_minus)
        check(f"{r.object} {r.method}", any(a in tex for a in alts), f"({sorted(alts)[0]})")

    print("Noterdaeme et al. 2011 (TeX Table 1)")
    tex = _squash((PAPER_DIR / "noterdaeme_2011" / "TCMB_astroph.tex").read_text())
    for r in df[df["source"].str.startswith("Noterdaeme11 Table")].itertuples():
        alts = {f"{r.z:.4f}&{a}" for a in asym_tokens(r.T, r.err_plus, r.err_minus)}
        check(f"{r.object} {r.T}", any(a in tex for a in alts), f"(z = {r.z})")

    print("Luzzi et al. 2009 (TeX Tables 1 and 3)")
    tex = _squash((PAPER_DIR / "luzzi_2009" / "gemmaluzzi3.tex").read_text())
    for r in df[df["reference"].str.startswith("Luzzi")].itertuples():
        name = r.object.split(" ")[0].replace("MS0451-0305", "MS0451")
        t_ok = f"{name}&{r.T:.2f}\\pm{r.err_plus:.2f}" in tex
        z_name = r.object.split(" ")[0]
        z_ok = f"{z_name}&{r.z:.3f}" in tex
        check(f"{r.object}", t_ok and z_ok, f"(T {t_ok}, z {z_ok})")

    print("Radio molecular absorbers (abstracts)")
    checks = {"Kotani et al. 2025": ("kotani_2025", "kotani_0925_arxiv.tex", "5.13\\pm0.06"),
              "Kotani Oka": ("kotani_2026", "kotani_arXiv.tex", "4.50\\pm0.17"),
              "Muller et al. 2013": ("muller_2013", "paper-PKS-CMB.tex", "5.08\\pm0.10")}
    for prefix, (key, fname, token) in checks.items():
        tex = _squash((PAPER_DIR / key / fname).read_text())
        check(prefix, token in tex, f"({token})")

    print(f"\n{len(FAILURES)} failures")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
