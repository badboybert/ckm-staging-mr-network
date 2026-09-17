# -*- coding: utf-8 -*-
"""Machine-check every factual assertion in RESPONSE_TO_INTERNAL_REVIEW_R5.md.

R5 responds to the round-5 FINAL PRE-SUBMISSION review (M1-M10 must-fix, R1-R10 recommended). Like
its predecessors, a response letter DESCRIBES work rather than being it, so it is not believed until
this script is green. Two inherited rules:

  * NOTHING THE PACKAGE CAN SUPPLY IS HARDCODED. Expected values are re-derived from the shipped
    package, from manifest/CANONICAL_FACTS.json, or by RE-RUNNING the gate the letter cites. Where the
    letter states a count (placeholders, R10 copy-edit occurrences, the R1-R5 / R6-R10 tallies) the
    count is PARSED from the letter and compared to the package, so "does the letter state what the
    package contains" cannot rot.
  * THE PACKAGE VERSION IS DERIVED from build_submission_v1.py — a verifier pointed at a package you
    did not just cut proves nothing.

Absence claims (M1/M2/R2-R5) are checked against the PACKAGE; the letter legitimately QUOTES the
phrases it removed, so those quotes must not be mistaken for the defect (we test SURF, never LET).
"Not applied" claims (R6-R10) are verified as still-present/absent so the letter's honesty is checked
both ways. Every check below is designed to be able to FAIL; add a mutation test before trusting it.

Run:  PYTHONIOENCODING=utf-8 python scripts/verify_rebuttal_r5.py
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

import io, os, re, sys, json, glob, subprocess

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = P4_ROOT
BASE = os.path.join(ROOT, "independent_build")
FIG = os.path.join(BASE, "figures")
SCR = os.path.join(BASE, "scripts")

_bsv = io.open(os.path.join(SCR, "build_submission_v1.py"), encoding="utf-8").read()
_m = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _m, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _m.group(1)
# LETTER and PKG are overridable so the mutation ladder (verify_rebuttal_r5_mutations.py) can point
# this verifier at a MUTATED COPY without touching the real letter or package. PKG_VERSION stays
# derived from the real builder, so the §1 version check still compares the letter to the true cut.
PKG = os.environ.get("CKM_R5_PKG") or os.path.join(ROOT, f"submission package {PKG_VERSION}")
LETTER = os.environ.get("CKM_R5_LETTER") or os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R5.md")

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
    except Exception as exc:                       # a reader that dies must be loud, never skipped
        print(f"  READER FAILED on {r}: {exc}")
        raise
ALLTXT = "\n".join(SURF.values())
F = json.load(open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))

# Declarations-only text (author-supplied placeholders live here + in the manuscript twin).
DECL_SCOPE = "\n".join(v for k, v in SURF.items() if k.startswith("02_cover_declarations/"))

OK, BAD = [], []


def check(label, cond, detail=""):
    (OK if cond else BAD).append(label)
    print(("  ok   " if cond else "  ✗ FAILED: ") + label + ("" if cond else f"  [{detail}]"))


def run(script):
    r = subprocess.run([sys.executable, os.path.join(SCR, script)], capture_output=True,
                       text=True, encoding="utf-8", errors="replace")
    return (r.stdout or "") + (r.stderr or "")


print(f"verify_rebuttal_r5 — letter vs {os.path.basename(PKG)}")
print(f"  {len(SURF)} shipped surfaces read "
      f"({sum(1 for k in SURF if k.endswith('.pdf'))} figure PDFs, "
      f"{sum(1 for k in SURF if k.endswith('.docx'))} docx, "
      f"{sum(1 for k in SURF if k.endswith('.xlsx'))} xlsx)\n")

# ---- 0. liveness: prove the readers work before trusting a single absence claim -----------------
check("§0 the package was fully read", len(SURF) >= 60, f"{len(SURF)} text-bearing")
check("§0 the docx reader recovered real content", "Additional file" in SURF.get(
    "01_manuscript/CKM_Paper4_manuscript.docx", ""))
check("§0 the xlsx reader recovered real content", "Waist circumference" in SURF.get(
    "04_supplementary/Supplementary_Tables.xlsx", ""))
_pdfs = [len(v) for k, v in SURF.items() if k.endswith(".pdf")]
check("§0 every figure PDF is in the plausible 200–6,000 char band",
      _pdfs and min(_pdfs) > 200 and max(_pdfs) < 6000, f"{min(_pdfs)}–{max(_pdfs)}")

# ---- 1. package identity, parsed from the letter -----------------------------------------------
_v = re.search(r"submission package (v[\d.]+)", LET)
check("§1 the letter names the package the build cuts", _v and _v.group(1) == PKG_VERSION,
      f"letter={_v and _v.group(1)} build={PKG_VERSION}")
_nf = re.search(r"submission package v[\d.]+` \((\d+) files\)", LET)
_actual = sum(1 for p in glob.glob(os.path.join(PKG, "**", "*"), recursive=True) if os.path.isfile(p))
check("§1 the file count in the letter matches the package",
      _nf and int(_nf.group(1)) == _actual, f"letter={_nf and _nf.group(1)} package={_actual}")

# ---- 2. MUST-FIX items M1-M10 -------------------------------------------------------------------
# M1: calibrated language — the superseded phrases are gone from the SHIPPED package.
for phr in ["largely holds", "does not manufacture", "directed causal graph"]:
    check(f"§2 M1: '{phr}' absent from every shipped file", phr not in ALLTXT)

# M2: Supplementary Figure S2 — no stale main-figure ref, no categorical annotation.
_s2 = SURF["04_supplementary/SupplFig2.pdf"]
check("§2 M2: no 'Fig 6' cross-reference on Supplementary Figure S2",
      not re.search(r"\bFig\.?\s*6\b", _s2))
check("§2 M2: no 'Non-causal' categorical annotation on S2", "Non-causal" not in _s2)

# M3: S6d legend no longer clipped — the full 'Bonferroni survives' string is present.
check("§2 M3: 'Bonferroni survives' is present in full on Supplementary Figure S6",
      "Bonferroni survives" in SURF["04_supplementary/SupplFig6.pdf"])

# M4: workbook title rows expanded (>30pt, wrapped) on the flagged sheets, no revision-history wording.
import openpyxl                                                            # noqa: E402
_wb = openpyxl.load_workbook(os.path.join(PKG, "04_supplementary", "Supplementary_Tables.xlsx"))
for sh in ["S14_overlap", "S17_HF_subtypes", "S18_WHRadjBMI", "S21_SBP_strict_primary"]:
    _h = _wb[sh].row_dimensions[1].height
    check(f"§2 M4: {sh} title row is expanded above the clipping 30pt", _h and _h > 30, f"height={_h}")
for phr in ["earlier implementation", "now the PRIMARY", "primary analysis replaces",
            "replace an earlier"]:
    check(f"§2 M4: revision-history phrase '{phr}' absent from the workbook",
          phr not in SURF["04_supplementary/Supplementary_Tables.xlsx"])

# M5: one physical file per declared Additional file — Source_Data.zip present; 4 declared.
check("§2 M5: Source_Data.zip is present in the package",
      os.path.exists(os.path.join(PKG, "04_supplementary", "Source_Data.zip")))
_ninv = len(re.findall(r"Additional file \d+\s*\|", SURF["02_cover_declarations/DECLARATIONS.md"]))
check("§2 M5: exactly four Additional files are declared", _ninv == 4, f"{_ninv} declared")

# M6: Additional-files table cannot split — every row cantSplit, header repeats.
import docx                                                                # noqa: E402
_dd = docx.Document(os.path.join(PKG, "01_manuscript", "CKM_Paper4_manuscript.docx"))
_af = next(t for t in _dd.tables
           if "Additional file" in " ".join(c.text for r in t.rows for c in r.cells))
_cant = all("cantSplit" in r._tr.xml for r in _af.rows)
_hdr = "tblHeader" in _af.rows[0]._tr.xml
check("§2 M6: every Additional-files row is set not to break across pages", _cant)
check("§2 M6: the Additional-files header row repeats", _hdr)

# M7: title page carries no internal production metrics.
_tp = SURF["02_cover_declarations/TITLE_PAGE.docx"]
for phr in ["Manuscript metrics", "short variant", "see ABSTRACT.md", "Main text total"]:
    check(f"§2 M7: '{phr}' absent from the journal-facing title page", phr not in _tp)

# M8: authors / ORCIDs / funding / competing interests filled; repo+DOI still pending.
_decl = SURF["02_cover_declarations/DECLARATIONS.md"]
for nm, orc in [("Son Tung Nguyen", "0009-0009-7546-7143"),
                ("Trang Huyen Nguyen", "0000-0002-4079-7486"),
                ("Ko-Ting Chen", "0000-0002-6795-5974"),
                ("Bertrand Chin-Ming Tan", "0000-0002-2218-7115")]:
    check(f"§2 M8: {nm} + ORCID {orc} on the title page", nm in _tp and orc in _tp)
_credit = _decl
_first_roles = "Formal analysis, Investigation, Data curation, Software, Visualization"
check("§2 M8: authors 1 and 2 (S.T.N., T.H.N.) share the first-author CRediT roles",
      "S.T.N." in _credit and "T.H.N." in _credit and _credit.count(_first_roles) >= 2)
check("§2 M8: the funders are named (NSTC + Chang Gung)",
      "National Science and Technology Council" in _decl and "Chang Gung Memorial Hospital" in _decl)
check("§2 M8: the competing-interests statement is present",
      "no competing interests" in _decl)
# R6 (2026-08-20): the repo/DOI placeholders were replaced by a truthful interim statement (finalised
# before publication); no bracketed placeholder survives in the declarations.
check("§2 M8: repository/DOI carried as an interim statement, no bracketed placeholder (R6 superseded R5)",
      "[GitHub URL" not in DECL_SCOPE
      and "provided prior to publication" in DECL_SCOPE
      and not re.search(r"\[AUTHOR-SUPPLIED", DECL_SCOPE))

# M9: README/AUTHOR_TODO agree; README does not list journal/citation/abstract as outstanding.
_rdme = io.open(os.path.join(PKG, "README_SUBMISSION.md"), encoding="utf-8").read()
check("§2 M9: README names Cardiovascular Diabetology as selected", "Cardiovascular Diabetology" in _rdme)
check("§2 M9: README does not present journal choice as outstanding",
      not re.search(r"journal choice[^.]*outstanding|outstanding[^.]*journal choice", _rdme, re.I))

# ---- 3. RECOMMENDED items R1-R10 ----------------------------------------------------------------
_abs = SURF["01_manuscript/ABSTRACT.md"]
# R1 applied — the two rewordings are in the abstract.
check("§3 R1: the reworded staging phrase is in the abstract",
      # round 7 re-led the Results paragraph with adiposity and tightened this clause
      "25 of 27 cross-stage edges ran from lower to higher stages" in _abs)
check("§3 R1: the reworded MVMR phrase is in the abstract",
      "a 1-SD higher body mass index (BMI) retained associations" in _abs)
# R2-R5 applied — each targeted phrase is gone from the package.
check("§3 R2: 'directed causal graph' is gone (checked under M1 too)", "directed causal graph" not in ALLTXT)
check("§3 R3: the categorical 'does not manufacture' claim is gone", "does not manufacture" not in ALLTXT)
check("§3 R4: 'Adiposity therefore acts' is gone", "Adiposity therefore acts" not in ALLTXT)
check("§3 R5: the failed-package troubleshooting sentence is gone",
      "could not be installed" not in ALLTXT and "mr.raps" not in ALLTXT)
# R6/R7 NOT applied — the letter says so, and honesty is checked both ways.
check("§3 R6: the on-panel 'native units' note is indeed absent from Figure 1 (not applied)",
      "not comparable" not in SURF["03_main_figures/Figure1.pdf"])
# R6 (2026-08-20): R7 applied — 'qhet' relabelled 'Q-minimisation' on Figure 2 panel d to match the prose.
check("§3 R7: 'qhet' relabelled 'Q-minimisation' on Figure 2 panel d (R6 superseded R5's deferral)",
      "qhet" not in SURF["03_main_figures/Figure2.pdf"].lower()
      and "q-minimisation" in SURF["03_main_figures/Figure2.pdf"].lower())
# R8 NOT applied — 'correlated-marker negative control' still used. Asserted on the Results section
# specifically (where it is used substantively): a presence check over ALLTXT is weaker — it passes on
# a stray copy in any surface, and cannot be mutation-tested by removing it from the manuscript.
check("§3 R8: 'correlated-marker negative control' is still present (not applied)",
      "correlated-marker negative control" in SURF["01_manuscript/RESULTS.md"])
# R10 NOT applied — the letter states the surviving occurrence counts; re-derive them.
_man_md = "\n".join(v for k, v in SURF.items()
                    if k.startswith("01_manuscript/") and k.endswith(".md"))
# R6 (2026-08-20): R10 applied — 'body-mass index' -> 'body mass index' and 'under-powered' ->
# 'underpowered' harmonised throughout; neither American form should remain in the manuscript.
for phr in ["body-mass index", "under-powered"]:
    _got = _man_md.count(phr)
    check(f"§3 R10: '{phr}' harmonised out of the manuscript sections (R6 superseded R5's deferral)",
          _got == 0, f"package={_got}")

# ---- 4. Taiwan Biobank consistency + data-access wording (§4) -----------------------------------
for k, lbl in [("01_manuscript/METHODS.md", "Methods"),
               ("02_cover_declarations/DECLARATIONS.md", "Declarations Availability"),
               ("04_supplementary/Supplementary_Tables.xlsx", "Supplementary Table S1"),
               ("04_supplementary/SUPPLEMENTARY_INFORMATION.md", "Supplementary Methods")]:
    check(f"§4 Taiwan Biobank is named in the {lbl}", "Taiwan Biobank" in SURF[k])
# The workbook reader joins each CELL on its own line, so a same-line regex over it cannot see a row;
# scan the S1 sheet's rows directly for one that carries both the trait and the cohort.
_s1 = _wb["S1_data_sources"]
_twb_row = any(
    any(isinstance(c, str) and "Waist circumference" in c for c in row) and
    any(isinstance(c, str) and "Taiwan Biobank" in c for c in row)
    for row in _s1.iter_rows(values_only=True))
check("§4 Supplementary Table S1 carries a Taiwan Biobank waist source row", _twb_row)
check("§4 the East-Asian data-access wording is calibrated (no 'no-DUA' claim)",
      "no institutional data-use agreement" not in ALLTXT
      and "European summary statistics are publicly available" in _decl)

# ---- 5. LLM disclosure — statement present, confirmation gate open (§5) --------------------------
_meth = SURF["01_manuscript/METHODS.md"]
check("§5 the Methods LLM statement names Anthropic Claude", "Anthropic Claude" in _meth)
# R6 (2026-08-20): the author-facing "[AUTHOR-SUPPLIED — confirm ...]" instruction was removed and the
# disclosure made publication-facing (Claude-only, author-confirmed).
check("§5 the LLM disclosure is publication-facing (R6: confirmation instruction removed)",
      "confirm that this statement covers every tool" not in _meth
      and "No AI system accessed raw data" in _meth)

# ---- 6. the summary tallies and placeholder count, parsed from the letter -----------------------
_m8 = re.search(r"Must-fix \(M1[–-]M10\) closed \| (\d+) of 10", LET)
check("§1 the letter's must-fix tally (8 of 10) is stated", _m8 and int(_m8.group(1)) == 8,
      f"letter={_m8 and _m8.group(1)}")
_r15 = re.search(r"\(R1[–-]R5\) applied \| (\d+) of 5", LET)
check("§1 the letter's R1-R5 tally (5 of 5) is stated", _r15 and int(_r15.group(1)) == 5,
      f"letter={_r15 and _r15.group(1)}")
_r610 = re.search(r"\(R6[–-]R10\) applied \| (\d+) of 5", LET)
check("§1 the letter's R6-R10 tally (0 of 5) is stated", _r610 and int(_r610.group(1)) == 0,
      f"letter={_r610 and _r610.group(1)}")
# The three remaining author-supplied field TYPES, derived from the package, must equal the letter's 3.
_types = set()
if "[GitHub URL" in DECL_SCOPE:
    _types.add("repo")
if re.search(r"Zenodo \(DOI \[AUTHOR-SUPPLIED\]\)", DECL_SCOPE):
    _types.add("doi")
if "confirm that this statement covers every tool" in ALLTXT:
    _types.add("llm")
# R6 (2026-08-20): all three tracked placeholders (repo URL, Zenodo DOI, LLM confirmation) were resolved
# — repo/DOI to a truthful interim statement, the LLM disclosure finalised — so none remain in the package.
check("§6 no tracked author-supplied placeholder remains in the package (R6 resolved all three)",
      len(_types) == 0, f"package={sorted(_types)}")
# No author/funding/COI/data-access placeholder must survive (those are claimed filled).
for phr in ["co-author", "affiliation", "CRediT taxonomy", "competing interests, unless",
            "grant numbers and funders", "institutional data-use agreements"]:
    check(f"§6 no surviving '{phr}' placeholder (the letter says these are filled)",
          not re.search(r"\[AUTHOR-SUPPLIED[^\]]*" + re.escape(phr), ALLTXT, re.I))

# ---- 7. §6 — the two gates the letter makes a countable claim about are RE-RUN, not quoted -------
# The mutation ladder sets CKM_R5_SKIP_RERUN: these two re-run the REAL package (their paths are
# absolute), so they neither see the temp mutation nor need to, and skipping them keeps each of the
# ladder's ~34 verifier invocations fast.
if not os.environ.get("CKM_R5_SKIP_RERUN"):
    _qa = run("qa_manuscript.py")
    _mqa = re.search(r"QA RESULTS: (\d+) PASS, (\d+) FAIL", _qa)
    check("§6 qa_manuscript.py returns 0 FAIL, as the letter's table states",
          _mqa and _mqa.group(2) == "0", f"actual={_mqa and _mqa.groups()}")
    _sc = run("scan_claims.py")
    check("§6 scan_claims.py is clean, as the letter's table states", "CLAIM SCAN CLEAN" in _sc)

# ---- 8. the letter must not over-claim: no minted DOI, no build-tree path in reader text ---------
check("§8 the letter does not claim a DOI is already minted",
      not re.search(r"10\.\d{4,}/zenodo", LET))
check("§8 the letter states the repository/DOI are pending",
      "inserted on issuance" in LET)

print("\n" + "=" * 96)
print(f"RESULT: {len(OK)} verified, {len(BAD)} FAILED")
if BAD:
    print("\nFAILED CLAIMS (the letter says something the package does not support):")
    for b in BAD:
        print("  x", b)
print("=" * 96)
sys.exit(1 if BAD else 0)
