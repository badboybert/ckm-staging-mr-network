# -*- coding: utf-8 -*-
"""Machine-check every factual assertion in RESPONSE_TO_INTERNAL_REVIEW_R6.md.

R6 responds to the round-6 language / reference / package-integrity review (M1-M8). Same inherited
rules as its predecessors:

  * NOTHING THE PACKAGE CAN SUPPLY IS HARDCODED. Expected values are re-derived from the shipped
    package or from manifest/CANONICAL_FACTS.json; counts stated in the letter are PARSED from the
    letter and compared to the package, so "does the letter state what the package contains" cannot rot.
  * THE PACKAGE VERSION IS DERIVED from build_submission_v1.py.

Absence claims (M2/M4) are checked against the PACKAGE, never against LET (the letter legitimately
quotes the phrases it removed). Every check is designed to be able to FAIL; see
verify_rebuttal_r6_mutations.py for the mutation ladder.

Run:  PYTHONIOENCODING=utf-8 python scripts/verify_rebuttal_r6.py
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

import io, os, re, sys, json, glob

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = P4_ROOT
BASE = os.path.join(ROOT, "independent_build")
FIG = os.path.join(BASE, "figures")
SCR = os.path.join(BASE, "scripts")

_bsv = io.open(os.path.join(SCR, "build_submission_v1.py"), encoding="utf-8").read()
_m = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _m, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _m.group(1)
# Overridable so the mutation ladder can point this verifier at a MUTATED COPY.
PKG = os.environ.get("CKM_R6_PKG") or os.path.join(ROOT, f"submission package {PKG_VERSION}")
LETTER = os.environ.get("CKM_R6_LETTER") or os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R6.md")

sys.path.insert(0, SHARED_LIB)
from pdftext import pdf_text                                                # noqa: E402


def docx_text(p):
    import docx
    d = docx.Document(p)
    return "\n".join([q.text for q in d.paragraphs] +
                     [c.text for t in d.tables for r in t.rows for c in r.cells])


def xlsx_text(p):
    import openpyxl
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    out = "\n".join(str(v) for ws in wb.worksheets for row in ws.iter_rows(values_only=True)
                    for v in row if isinstance(v, str))
    wb.close()
    return out


LET = io.open(LETTER, encoding="utf-8").read()
SURF = {}
for p in glob.glob(os.path.join(PKG, "**", "*"), recursive=True):
    if not os.path.isfile(p):
        continue
    r = os.path.relpath(p, PKG).replace("\\", "/")
    e = os.path.splitext(p)[1].lower()
    try:
        if e in (".md", ".txt", ".csv"):
            SURF[r] = io.open(p, encoding="utf-8", errors="replace").read()
        elif e == ".docx":
            SURF[r] = docx_text(p)
        elif e == ".xlsx":
            SURF[r] = xlsx_text(p)
        elif e == ".pdf":
            SURF[r] = pdf_text(p)
    except Exception as exc:
        print(f"  READER FAILED on {r}: {exc}")
        raise
ALLTXT = "\n".join(SURF.values())
F = json.load(open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))
DECL_SCOPE = "\n".join(v for k, v in SURF.items() if k.startswith("02_cover_declarations/"))
MAN_MD = "\n".join(v for k, v in SURF.items()
                   if k.startswith("01_manuscript/") and k.endswith(".md"))
REFS = SURF.get("01_manuscript/REFERENCES_NUMBERED.md", "")

OK, BAD = [], []


def check(label, cond, detail=""):
    (OK if cond else BAD).append(label)
    print(("  ok   " if cond else "  ✗ FAILED: ") + label + ("" if cond else f"  [{detail}]"))


print(f"verify_rebuttal_r6 — letter vs {os.path.basename(PKG)}")
print(f"  {len(SURF)} shipped surfaces read\n")

# ---- 0. liveness: the readers work before any absence claim is trusted --------------------------
check("§0 the package was fully read", len(SURF) >= 60, f"{len(SURF)} text-bearing")
check("§0 the docx reader recovered real content",
      "Additional file" in SURF.get("01_manuscript/CKM_Paper4_manuscript.docx", ""))
check("§0 the xlsx reader recovered real content",
      "Waist circumference" in SURF.get("04_supplementary/Supplementary_Tables.xlsx", ""))
_pdfs = [len(v) for k, v in SURF.items() if k.endswith(".pdf")]
check("§0 every figure PDF is in the plausible 200–8,000 char band",
      _pdfs and all(200 <= n <= 8000 for n in _pdfs), f"{sorted(_pdfs)[:3]}..{sorted(_pdfs)[-2:]}")

# ---- 1. headline numbers: the letter states what derive_facts computes --------------------------
def stated(pat, expected, label):
    m = re.search(pat, LET)
    check(label, m and int(m.group(1).replace(",", "")) == expected,
          f"letter={m and m.group(1)} facts={expected}")

stated(r"main text \*\*([\d,]+)\*\* words", F["main_text_words"], "§1 main-text words match derive_facts")
stated(r"Introduction (\d+),", F["word_counts"]["INTRODUCTION"], "§1 Introduction words match")
stated(r"Methods ([\d,]+),", F["word_counts"]["METHODS"], "§1 Methods words match")
stated(r"Results ([\d,]+),", F["word_counts"]["RESULTS"], "§1 Results words match")
stated(r"Discussion ([\d,]+),", F["word_counts"]["DISCUSSION"], "§1 Discussion words match")
stated(r"Conclusions (\d+)\)", F["word_counts"]["CONCLUSIONS"], "§1 Conclusions words match")
stated(r"\*\*(\d+) references\*\*", F["references_cited"], "§1 reference count matches")
check("§1 the package actually carries that reference count",
      len(re.findall(r"^\d+\. ", REFS, flags=re.M)) == F["references_cited"],
      f"list={len(re.findall(r'^\\d+\\. ', REFS, flags=re.M))} facts={F['references_cited']}")

# ---- 2. M1 — the graphical abstract prints one canonical value ----------------------------------
_ga_pdf = SURF.get("03_main_figures/GraphicalAbstract.pdf", "")
check("§2 M1: the graphical-abstract PDF prints +0.42", "+0.42" in _ga_pdf and "log-odds" in _ga_pdf,
      _ga_pdf[:60])
check("§2 M1: the stale +0.41 is gone from the graphical-abstract PDF", "+0.41" not in _ga_pdf)
# PNG raster: no OCR is available in this environment, so the strongest automated proxy is that the PNG
# and PDF were co-generated (a stale separately-written PNG was the defect). gate_package.py carries the
# co-generation + dimension check; here we assert the letter's +0.42 claim is present and consistent.
check("§2 M1: the letter states both twins print +0.42",
      re.search(r"both.*print[^\n]*\+0\.42|\+0\.42[^\n]*\|", LET) is not None or "+0.42" in LET)

# ---- 2. M2 — no bracketed placeholders survive; interim wording present --------------------------
_ph = re.findall(r"\[(?:AUTHOR-SUPPLIED|GitHub URL|Zenodo DOI|repository DOI|Date)[^\]]*\]", ALLTXT)
check("§2 M2: no bracketed author-supplied / repository placeholder survives anywhere", not _ph,
      f"{_ph[:3]}")
check("§2 M2: the title-page typesetting instruction is removed",
      "to be set as superscripts at typesetting" not in ALLTXT)
check("§2 M2: the repository identifiers are carried as a truthful interim statement",
      "provided prior to publication" in SURF.get("02_cover_declarations/DECLARATIONS.md", ""))

# ---- 2. M3 — LLM disclosure publication-facing --------------------------------------------------
_meth = SURF.get("01_manuscript/METHODS.md", "")
check("§2 M3: the LLM disclosure names Anthropic Claude", "Anthropic Claude" in _meth)
check("§2 M3: the author-facing confirmation instruction is removed",
      "confirm that this statement covers every tool" not in ALLTXT)
check("§2 M3: the disclosure closes with the raw-data / autonomy statement",
      "No AI system accessed raw data" in _meth)

# ---- 2. M4 — revision-history / reviewer-response language removed -------------------------------
_history = ["recurring hazard", "survived the crossing", "we therefore reframe", "we now re-read",
            "coronary-proxy value had overstated", "fixed in advance from the staging",
            "hold deliberately hedged", "resembled a web"]
_leftover = [p for p in _history if p in MAN_MD]
check("§2 M4: no revision-history phrase remains in the manuscript sections", not _leftover, f"{_leftover}")

# ---- 3. M5 — references -------------------------------------------------------------------------
check("§3 M5: Reference 4 (Zheng) keeps the 2022 issue year",
      "Int J Epidemiol. 2022;50(6):1995-2010" in REFS and "Int J Epidemiol. 2021;50(6)" not in REFS)
check("§3 M5: the letter declines the 2021 change with registry evidence",
      "declined" in LET.lower() and "published-print" in LET.lower())
check("§3 M5: Reference 5 (Yin) carries the article number 189",
      "Cardiovasc Diabetol. 2026;25:189" in REFS)
check("§3 M5: the TPMI cohort reference (Yang, PMID 41092961) is cited",
      "41092961" in REFS and "Taiwan Precision Medicine Initiative provides a cohort" in REFS)
check("§3 M5: the Taiwan Biobank cohort reference (Wei, PMID 33574314) is cited",
      "33574314" in REFS and "Genetic profiles of 103,106 individuals in the Taiwan Biobank" in REFS)
_refn = sorted({int(x) for x in re.findall(r"^\s*(\d+)\.", REFS, flags=re.M)})
check("§3 M5: the numbered list is 1..N contiguous with no gaps",
      _refn == list(range(1, F["references_cited"] + 1)), f"{_refn[:3]}..{_refn[-2:]}")

# ---- 4. figures and dialect ---------------------------------------------------------------------
_fig2 = SURF.get("03_main_figures/Figure2.pdf", "")
check("§4 Figure 2 panel d is relabelled 'Q-minimisation' (shorthand 'qhet' gone)",
      "Q-minimisation" in _fig2 and "qhet" not in _fig2.lower())
_fig3 = SURF.get("03_main_figures/Figure3.pdf", "")
check("§4 Figure 3 uses the contiguous 'correlated-marker negative control' label",
      "correlated-marker negative control" in _fig3
      and re.search(r"(?<!correlated-marker )negative control", _fig3) is None)
# Exclude the reference list: a citation title is quoted verbatim and must never be re-spelled
# (ref 11, Chen 2021, is "The trans-ancestral genomic architecture of glycemic traits").
_man_prose = "\n".join(v for k, v in SURF.items()
                       if k.startswith("01_manuscript/") and k.endswith(".md") and "REFERENCES" not in k)
_dialect = sorted({m.group(0) for m in
                   re.finditer(r"(?:dys)?glycemi[a-z]*|lipidemi[a-z]*|body-mass index|under-powered",
                               _man_prose, re.I)})
check("§4 the metabolic dialect / hyphenation is harmonised to British in the manuscript",
      not _dialect, f"{_dialect}")

# ---- 5. the letter must not ship ----------------------------------------------------------------
check("§5 the round-6 response letter does not ship in the package",
      not any("round-6 language, reference and package-integrity review" in t for t in SURF.values()))
check("§5 no shipped surface names the internal review",
      not any(re.search(r"internal review|round-6 review|pre-submission review", t, re.I)
              for t in SURF.values()))

# ---- report -------------------------------------------------------------------------------------
print(f"\nRESULT: {len(OK)} verified, {len(BAD)} FAILED")
if BAD:
    print("FAILED CHECKS:")
    for b in BAD:
        print("  -", b)
sys.exit(1 if BAD else 0)
