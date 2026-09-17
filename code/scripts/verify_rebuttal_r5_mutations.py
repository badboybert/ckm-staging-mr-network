# -*- coding: utf-8 -*-
"""Mutation-test ladder for verify_rebuttal_r5.py.

A check that stays green under a defect that violates it is vacuous. For every distinct check the
verifier makes, this ladder injects one targeted defect into a COPY of the letter or package that
SHOULD trip that check, runs verify_rebuttal_r5.py against the copy, and asserts THE NAMED CHECK goes
red. It then restores and moves on. Two disciplines it enforces on itself:

  * EVERY MUTATION MUST ACTUALLY CHANGE THE INPUT. A regex/edit that is a no-op (its target string
    has moved or been renamed) is refused loudly — otherwise the ladder reports a MISS that looks like
    a broken gate when the test itself is stale.
  * THE BASELINE COPY MUST BE GREEN FIRST. If the untouched temp copy is not 0-failed, the workspace
    is wrong and no mutation result can be trusted.

The verifier is pointed at the copy through CKM_R5_LETTER / CKM_R5_PKG, and CKM_R5_SKIP_RERUN makes
each of the ~30 invocations fast (qa/scan_claims re-run the real package and are irrelevant here).

Run:  PYTHONIOENCODING=utf-8 python scripts/verify_rebuttal_r5_mutations.py
"""
# --- portable roots -----------------------------------------------------------------------------
# Injected by scripts/build_repo.py. The working tree hardcoded an absolute local path; the deposit
# resolves it from CKM_P4_BASE, or from this file's own location (repo/code/ plays the role the
# working tree called independent_build/). Raw GWAS summary statistics are NOT deposited: set
# CKM_ROOT to wherever you obtained them if you intend to re-run the upstream extraction steps.
import os as _os
P4_BASE = _os.environ.get("CKM_P4_BASE") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
P4_ROOT = _os.environ.get("CKM_P4_ROOT") or _os.path.dirname(P4_BASE)
CKM_ROOT = _os.environ.get("CKM_ROOT") or _os.path.dirname(P4_ROOT)
SHARED_LIB = _os.environ.get("CKM_SHARED_LIB") or _os.path.join(CKM_ROOT, "_shared")
# ------------------------------------------------------------------------------------------------

import io, os, re, sys, shutil, tempfile, subprocess

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, SHARED_LIB)
from pdftext import pdf_text                                                # noqa: E402

ROOT = P4_ROOT
BASE = os.path.join(ROOT, "independent_build")
SCR = os.path.join(BASE, "scripts")
_bsv = io.open(os.path.join(SCR, "build_submission_v1.py"), encoding="utf-8").read()
PKG_VERSION = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M).group(1)
R_PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")
R_LETTER = os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R5.md")

WORK = tempfile.mkdtemp(prefix="r5mut_")
TMP_PKG = os.path.join(WORK, "pkg")
TMP_LETTER = os.path.join(WORK, "letter.md")
shutil.copytree(R_PKG, TMP_PKG)
shutil.copy(R_LETTER, TMP_LETTER)


def restore(rel):
    if rel == "letter":
        shutil.copy(R_LETTER, TMP_LETTER)
    else:
        shutil.copy(os.path.join(R_PKG, rel), os.path.join(TMP_PKG, rel))


# ---------- mutation primitives (each REFUSES a no-op) -------------------------------------------
def _sub(path, find, repl, count=0):
    t = io.open(path, encoding="utf-8").read()
    t2 = re.sub(find, repl, t, count=count)
    assert t2 != t, f"NO-OP mutation (target gone?): {find!r} in {os.path.basename(path)}"
    io.open(path, "w", encoding="utf-8").write(t2)


def L(find, repl, count=0):          # letter
    _sub(TMP_LETTER, find, repl, count)


def M(rel, find, repl, count=0):     # a text surface (.md) in the package
    _sub(os.path.join(TMP_PKG, rel), find, repl, count)


def M_append(rel, text):             # add text to a .md — trips an ABSENCE-in-ALLTXT check
    io.open(os.path.join(TMP_PKG, rel), "a", encoding="utf-8").write(text)


def M_del(rel):
    p = os.path.join(TMP_PKG, rel)
    assert os.path.exists(p), f"NO-OP delete (already gone): {rel}"
    os.remove(p)


def X(fn):                           # mutate the supplementary workbook
    import openpyxl
    p = os.path.join(TMP_PKG, "04_supplementary", "Supplementary_Tables.xlsx")
    wb = openpyxl.load_workbook(p)
    fn(wb)
    wb.save(p)


def D(rel, fn):                      # mutate a .docx
    import docx
    p = os.path.join(TMP_PKG, rel)
    d = docx.Document(p)
    fn(d)
    d.save(p)


def PDF_add(rel, text):              # insert visible text into a figure PDF (trips an ABSENCE check)
    import fitz
    p = os.path.join(TMP_PKG, rel)
    d = fitz.open(p)
    d[0].insert_text(fitz.Point(72, 72), text, fontsize=11)
    d.save(p, incremental=True, encryption=fitz.PDF_ENCRYPT_KEEP)
    d.close()
    assert text.split()[0] in pdf_text(p), f"PDF add did not take in {rel}"


def PDF_swap(rel, src_rel):          # overwrite a PDF with another real one (trips a PRESENCE check)
    shutil.copy(os.path.join(TMP_PKG, src_rel), os.path.join(TMP_PKG, rel))


def docx_replace(d, find, repl):
    # Replace in EVERY matching paragraph, not just the first — a name or phrase recurs (the author
    # line AND the ORCID line; a term used in Results AND Methods), and a surviving copy is a MISS.
    paras = d.paragraphs + [c.paragraphs[0] for t in d.tables for r in t.rows for c in r.cells if c.paragraphs]
    hit = False
    for para in paras:
        if find in para.text:
            new = para.text.replace(find, repl)
            for r in para.runs:
                r.text = ""
            if para.runs:
                para.runs[0].text = new
            else:
                para.add_run(new)
            hit = True
    assert hit, f"NO-OP docx replace (not found): {find}"


def docx_uncantsplit(d):
    from docx.oxml.ns import qn
    hit = False
    for t in d.tables:
        if any("Additional file" in c.text for r in t.rows for c in r.cells):
            for r in t.rows:
                trPr = r._tr.find(qn("w:trPr"))
                if trPr is not None:
                    cs = trPr.find(qn("w:cantSplit"))
                    if cs is not None:
                        trPr.remove(cs)
                        hit = True
    assert hit, "NO-OP: no cantSplit found to remove"


# ---------- the ladder: (expected-failing check label, [files to restore], apply) ----------------
_ARROWS = r"body-mass index"
MUT = [
    # ---- letter-side ----
    ("§1 the letter names the package", ["letter"],
     lambda: L(r"submission package v[\d.]+", "submission package v9.9")),
    ("§1 the file count in the letter", ["letter"], lambda: L(r"\(99 files\)", "(98 files)", 1)),
    ("§1 the letter's must-fix tally", ["letter"], lambda: L(r"8 of 10", "7 of 10")),
    ("§1 the letter's R1-R5 tally", ["letter"], lambda: L(r"5 of 5", "4 of 5")),
    ("§1 the letter's R6-R10 tally", ["letter"], lambda: L(r"0 of 5", "1 of 5")),
    ("§3 R10: the letter's 'body-mass index' count", ["letter"],
     lambda: L(_ARROWS + r'(["”]?) \(×5\)', _ARROWS + r"\1 (×9)")),
    ("§3 R10: the letter's 'under-powered' count", ["letter"],
     lambda: L(r'under-powered(["”]?) \(×2\)', r"under-powered\1 (×9)")),
    ("§6 the letter's remaining-placeholder count", ["letter"], lambda: L(r"placeholders \| 3", "placeholders | 4")),
    ("§8 the letter does not claim a DOI is already minted", ["letter"],
     lambda: L(r"(## 1\. Summary)", r"Archived at 10.5281/zenodo.99999.\n\n\1", 1)),
    ("§8 the letter states the repository/DOI are pending", ["letter"],
     lambda: L(r"inserted on issuance", "inserted eventually")),
    # ---- package .md ----
    ("§2 M1: 'largely holds' absent", ["01_manuscript/RESULTS.md"],
     lambda: M_append("01_manuscript/RESULTS.md", "\nThe architecture largely holds.\n")),
    ("§2 M5: exactly four Additional files", ["02_cover_declarations/DECLARATIONS.md"],
     lambda: M("02_cover_declarations/DECLARATIONS.md", r"Additional file 4 \|", "Additional file X |", 1)),
    ("§2 M8: the funders are named", ["02_cover_declarations/DECLARATIONS.md"],
     lambda: M("02_cover_declarations/DECLARATIONS.md", "National Science and Technology Council", "Redacted Agency")),
    ("§2 M8: the competing-interests statement is present", ["02_cover_declarations/DECLARATIONS.md"],
     lambda: M("02_cover_declarations/DECLARATIONS.md", "no competing interests", "several competing interests")),
    ("§2 M9: README names Cardiovascular Diabetology", ["README_SUBMISSION.md"],
     lambda: M("README_SUBMISSION.md", "Cardiovascular Diabetology", "Some Other Journal")),
    ("§3 R1: the reworded staging phrase", ["01_manuscript/ABSTRACT.md"],
     lambda: M("01_manuscript/ABSTRACT.md", "25 of 27 cross-stage edges pointed from lower to higher stages", "the edges were concordant")),
    ("§3 R2: 'directed causal graph' is gone", ["01_manuscript/RESULTS.md"],
     lambda: M_append("01_manuscript/RESULTS.md", "\nWe built a directed causal graph.\n")),
    ("§3 R8: 'correlated-marker negative control' is still present", ["01_manuscript/RESULTS.md"],
     lambda: M("01_manuscript/RESULTS.md", "correlated-marker negative control", "calibration marker")),
    ("§4 Taiwan Biobank is named in the Methods", ["01_manuscript/METHODS.md"],
     lambda: M("01_manuscript/METHODS.md", "Taiwan Biobank", "a Taiwan cohort")),
    ("§4 the East-Asian data-access wording is calibrated", ["02_cover_declarations/DECLARATIONS.md"],
     lambda: M_append("02_cover_declarations/DECLARATIONS.md", "\nAccess required no institutional data-use agreement.\n")),
    ("§5 the Methods LLM statement names Anthropic Claude", ["01_manuscript/METHODS.md"],
     lambda: M("01_manuscript/METHODS.md", "Anthropic Claude", "an assistant")),
    ("§5 the LLM-disclosure confirmation gate is still open", ["01_manuscript/METHODS.md"],
     lambda: M("01_manuscript/METHODS.md", r"\[AUTHOR-SUPPLIED — confirm that this statement covers every tool[^\]]*\]", "")),
    ("§6 no surviving 'co-author' placeholder", ["02_cover_declarations/TITLE_PAGE.md"],
     lambda: M_append("02_cover_declarations/TITLE_PAGE.md", "\n[AUTHOR-SUPPLIED — co-author list to confirm]\n")),
    # ---- workbook (.xlsx) ----
    ("§2 M4: S14_overlap title row is expanded", ["04_supplementary/Supplementary_Tables.xlsx"],
     lambda: X(lambda wb: setattr(wb["S14_overlap"].row_dimensions[1], "height", 20))),
    ("§2 M4: revision-history phrase 'now the PRIMARY' absent", ["04_supplementary/Supplementary_Tables.xlsx"],
     lambda: X(lambda wb: wb["Contents"].append(["now the PRIMARY analysis"]))),
    ("§4 Supplementary Table S1 carries a Taiwan Biobank waist source row",
     ["04_supplementary/Supplementary_Tables.xlsx"],
     lambda: X(lambda wb: [setattr(c, "value", "redacted")
                           for row in wb["S1_data_sources"].iter_rows() for c in row
                           if isinstance(c.value, str) and "Taiwan Biobank" in c.value])),
    # ---- title page / manuscript (.docx) ----
    ("§2 M7: 'Manuscript metrics' absent from the journal-facing title page",
     ["02_cover_declarations/TITLE_PAGE.docx"],
     lambda: D("02_cover_declarations/TITLE_PAGE.docx", lambda d: d.add_paragraph("Manuscript metrics"))),
    ("§2 M8: Son Tung Nguyen + ORCID 0009-0009-7546-7143 on the title page",
     ["02_cover_declarations/TITLE_PAGE.docx"],
     lambda: D("02_cover_declarations/TITLE_PAGE.docx", lambda d: docx_replace(d, "Son Tung Nguyen", "Author Redacted"))),
    ("§2 M6: every Additional-files row is set not to break across pages",
     ["01_manuscript/CKM_Paper4_manuscript.docx"],
     lambda: D("01_manuscript/CKM_Paper4_manuscript.docx", docx_uncantsplit)),
    # ---- figure PDFs ----
    ("§2 M2: no 'Fig 6' cross-reference on Supplementary Figure S2", ["04_supplementary/SupplFig2.pdf"],
     lambda: PDF_add("04_supplementary/SupplFig2.pdf", "Non-causal Fig 6")),
    ("§2 M3: 'Bonferroni survives' is present in full", ["04_supplementary/SupplFig6.pdf"],
     lambda: PDF_swap("04_supplementary/SupplFig6.pdf", "04_supplementary/SupplFig1.pdf")),
    ("§3 R7: the shorthand 'qhet' indeed remains", ["03_main_figures/Figure2.pdf"],
     lambda: PDF_swap("03_main_figures/Figure2.pdf", "03_main_figures/Figure1.pdf")),
    ("§3 R6: the on-panel 'native units' note is indeed absent", ["03_main_figures/Figure1.pdf"],
     lambda: PDF_add("03_main_figures/Figure1.pdf", "magnitudes are not comparable across rows")),
]


def run_verifier():
    env = dict(os.environ, CKM_R5_LETTER=TMP_LETTER, CKM_R5_PKG=TMP_PKG, CKM_R5_SKIP_RERUN="1",
               PYTHONIOENCODING="utf-8")
    r = subprocess.run([sys.executable, os.path.join(SCR, "verify_rebuttal_r5.py")],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", env=env)
    return (r.stdout or "") + (r.stderr or "")


def failed_labels(out):
    return set(re.findall(r"✗ FAILED: (.+?)(?:  \[|$)", out, flags=re.M))


TRIPPED, MISSED = [], []
try:
    # baseline: the untouched copy must be fully green, or nothing below can be trusted
    base = run_verifier()
    _bl = re.search(r"RESULT: (\d+) verified, (\d+) FAILED", base)
    assert _bl and _bl.group(2) == "0", "BASELINE COPY IS NOT GREEN — workspace is wrong:\n" + base[-1500:]
    print(f"baseline: {_bl.group(1)} verified, 0 FAILED (clean copy)\n")

    for label, files, apply in MUT:
        for f in files:
            restore(f)
        apply()                       # raises on a no-op mutation
        out = run_verifier()
        fl = failed_labels(out)
        hit = any(label in x for x in fl)
        (TRIPPED if hit else MISSED).append(label)
        print(("  ✓ tripped  " if hit else "  ✗ MISSED  ") + label
              + ("" if hit else f"   [failing checks were: {sorted(fl) or 'NONE'}]"))
        for f in files:               # leave the workspace clean for the next mutation
            restore(f)
finally:
    shutil.rmtree(WORK, ignore_errors=True)

print("\n" + "=" * 96)
print(f"LADDER: {len(TRIPPED)} of {len(MUT)} mutations tripped their named check")
if MISSED:
    print("\nVACUOUS OR MIS-TARGETED CHECKS (a defect that violates them did NOT make them fail):")
    for m in MISSED:
        print("  x", m)
print("=" * 96)
sys.exit(1 if MISSED else 0)
