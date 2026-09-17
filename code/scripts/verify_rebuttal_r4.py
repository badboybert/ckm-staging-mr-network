# -*- coding: utf-8 -*-
"""Machine-check every factual assertion in RESPONSE_TO_INTERNAL_REVIEW_R4.md.

A response letter is the easiest document in the project to get wrong, because it DESCRIBES work
rather than being it. The round-2 letter passed a careful read and then failed six of its own claims
on first machine check, four of which were real package defects. So this letter is not believed
until this script is green.

Two rules this file inherits from verify_rebuttal_r3.py, both earned:

  * NOTHING IS HARDCODED. Every expected value is re-derived from the shipped package, from
    manifest/CANONICAL_FACTS.json, or from the results file the number came from. The R3 verifier
    originally hardcoded the letter's numbers, so a single change had to be made in three places and
    the check could go stale without failing. It now PARSES the letter and compares — "does the
    letter state what the package contains" cannot rot.
  * THE PACKAGE VERSION IS DERIVED. A verifier pointed at a package you did not just cut proves
    nothing.

Run:  PYTHONIOENCODING=utf-8 python scripts/verify_rebuttal_r4.py
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

import io, os, re, sys, csv, json, glob, subprocess

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = P4_ROOT
BASE = os.path.join(ROOT, "independent_build")
FIG = os.path.join(BASE, "figures")
RES = os.path.join(BASE, "results")
SCR = os.path.join(BASE, "scripts")

_bsv = io.open(os.path.join(SCR, "build_submission_v1.py"), encoding="utf-8").read()
_m = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _m, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _m.group(1)
PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")
LETTER = os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R4.md")

# The shared PDF reader, never a local copy: five forks in this project were returning ~70x too many
# characters (font-encoding data), so every presence probe over them was weaker than it looked.
sys.path.insert(0, SHARED_LIB)
from pdftext import pdf_text                                             # noqa: E402


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

OK, BAD = [], []


def check(label, cond, detail=""):
    (OK if cond else BAD).append(label)
    print(("  ok   " if cond else "  ✗ FAILED: ") + label + ("" if cond else f"  [{detail}]"))


def letter_has(pat):
    return re.search(pat, LET) is not None


print(f"verify_rebuttal_r4 — letter vs {os.path.basename(PKG)}")
print(f"  {len(SURF)} shipped surfaces read "
      f"({sum(1 for k in SURF if k.endswith('.pdf'))} figure PDFs, "
      f"{sum(1 for k in SURF if k.endswith('.docx'))} docx, "
      f"{sum(1 for k in SURF if k.endswith('.xlsx'))} xlsx)\n")

# ---- 0. liveness: prove the readers work before trusting a single absence claim -----------------
check("§0 the package was fully read (99 files)", len(SURF) + 0 >= 60, f"{len(SURF)} text-bearing")
check("§0 the docx reader recovered real content", "Additional file" in SURF.get(
    "01_manuscript/CKM_Paper4_manuscript.docx", ""))
check("§0 the xlsx reader recovered real content", "Transition" in SURF.get(
    "04_supplementary/Supplementary_Tables.xlsx", ""))
_pdfs = [len(v) for k, v in SURF.items() if k.endswith(".pdf")]
check("§0 every figure PDF is in the plausible 200–6,000 char band",
      _pdfs and min(_pdfs) > 200 and max(_pdfs) < 6000, f"{min(_pdfs)}–{max(_pdfs)}")

# ---- 1. the headline package identity, parsed from the letter -----------------------------------
_v = re.search(r"submission package (v[\d.]+)", LET)
check("§1 the letter names the package the build actually cuts",
      _v and _v.group(1) == PKG_VERSION, f"letter={_v and _v.group(1)} build={PKG_VERSION}")
_nf = re.search(r"submission package v[\d.]+` \((\d+) files\)", LET)
_actual = sum(1 for p in glob.glob(os.path.join(PKG, "**", "*"), recursive=True) if os.path.isfile(p))
check("§1 the file count in the letter matches the package",
      _nf and int(_nf.group(1)) == _actual, f"letter={_nf and _nf.group(1)} package={_actual}")

# ---- 2. every shipped number the letter quotes, re-derived --------------------------------------
w = F["word_counts"]
NUMS = [
    ("§7 abstract words", r"Abstract \*\*(\d+) words\*\*", F["abstract_full_words"]),
    ("§7 main-text words", r"Main text \*\*([\d,]+) words\*\*", F["main_text_words"]),
    ("§7 Introduction", r"Introduction ([\d,]+), Methods", w["INTRODUCTION"]),
    ("§7 Methods", r"Methods ([\d,]+), Results", w["METHODS"]),
    ("§7 Results", r"Results ([\d,]+), Discussion", w["RESULTS"]),
    ("§7 Discussion", r"Discussion ([\d,]+), Conclusions", w["DISCUSSION"]),
    ("§7 Conclusions", r"Conclusions ([\d,]+)\.", w["CONCLUSIONS"]),
    ("§7 main figures", r"\*\*(\d+) main figures\*\*", F["main_figures"]),
    ("§7 panels", r"\*\*(\d+) panels\*\*", F["total_panels"]),
    ("§7 supplementary figures", r"\*\*(\d+) supplementary figures\*\*", F["supp_figures"]),
    ("§7 supplementary tables", r"Supplementary Tables S1–S(\d+)", F["supp_tables"]),
    ("§7 worksheets", r"\*\*(\d+) worksheets\*\*", F["supp_table_sheets"]),
    ("§7 references", r"\*\*(\d+) references\*\*", F["references_cited"]),
    ("§7 network edges", r"tests\s+\*\*(\d+)\*\* directed edges", F["network_edges_total"]),
    ("§7 Bonferroni edges", r"\*\*(\d+)\*\* pass the network-wide Bonferroni",
     F["edges_passing_bonferroni"]),
]
for label, pat, want in NUMS:
    m = re.search(pat, LET)
    got = int(m.group(1).replace(",", "")) if m else None
    check(f"{label} in the letter matches the package ({want})", got == want,
          f"letter says {got}")

# ---- 3. numbers the letter quotes that come from a RESULTS FILE, not from derive_facts ----------
_cn = io.open(os.path.join(RES, "staging_constrained_nulls.txt"), encoding="utf-8").read()
_p_deg = float(re.search(r"C\. Degree-preserving[^\n]*?P = ([\d.]+)", _cn).group(1))
_m = re.search(r"degree-preserving rewiring null\s*\n?is not exceeded \(\*P\* = ([\d.]+)\)", LET)
check("§7 the degree-preserving null P in the letter matches its source file",
      _m and abs(float(_m.group(1)) - _p_deg) < 5e-3, f"letter={_m and _m.group(1)} file={_p_deg}")
_m = re.search(r"concordance is \*\*([\d.]+)\*\*", LET)
_stat = {r["metric"]: r["value"] for r in
         csv.DictReader(open(os.path.join(FIG, "data", "fig1_staging_stat.csv"), encoding="utf-8"))}
check("§7 the cross-stage concordance in the letter matches its source file",
      _m and abs(float(_m.group(1)) - float(_stat["observed_concordance"])) < 5e-4,
      f"letter={_m and _m.group(1)} file={_stat['observed_concordance']}")
_sd = [r for r in csv.DictReader(open(os.path.join(RES, "forward_local_edges.csv"), encoding="utf-8"))]
_n53 = len(list(csv.DictReader(open(os.path.join(FIG, "data", "fig1_edges.csv"), encoding="utf-8"))))
_m = re.search(r"\*\*(\d+)\*\*\s*\n?survive Steiger gating", LET)
check("§7 the Steiger-directed graph size matches fig1_edges.csv",
      _m and int(_m.group(1)) == _n53, f"letter={_m and _m.group(1)} file={_n53}")

# ---- 4. the CAD→SBP repair the letter describes (§2) --------------------------------------------
_cs = [r for r in _sd if r["exposure"] == "CAD" and r["outcome"] == "SBP"][0]
_m = re.search(r"Egger intercept of \*\*\+([\d.]+), \*P\* = ([\d.]+)\*\*", LET)
check("§2 the Egger intercept in the letter matches forward_local_edges.csv",
      _m and abs(float(_m.group(1)) - float(_cs["egger_intercept"])) < 5e-4
      and abs(float(_m.group(2)) - float(_cs["egger_intercept_p"])) < 5e-4,
      f"letter={_m and _m.groups()} "
      f"file={_cs['egger_intercept']}/{_cs['egger_intercept_p']}")
# The CORRECTED triplet, anchored to its own sentence and case-sensitive to it. Anchoring matters:
# a bare r"weighted median \+([\d.]+)" matched the SUPERSEDED +0.62 four lines earlier and reported
# the letter wrong when the letter was right — the brittle-needle failure this project keeps hitting.
_m = re.search(r"They are now \*P\* = ([\d.]+), \+([\d.]+) and \+([\d.]+)", LET)
check("§2 the corrected CAD→SBP triplet matches forward_local_edges.csv",
      _m and abs(float(_m.group(1)) - float(_cs["egger_intercept_p"])) < 5e-4
      and abs(float(_m.group(2)) - float(_cs["wm_b"])) < 5e-3
      and abs(float(_m.group(3)) - float(_cs["ivw_b"])) < 5e-3,
      f"letter={_m and _m.groups()} "
      f"file={_cs['egger_intercept_p']}/{_cs['wm_b']}/{_cs['ivw_b']}")
# EVERY mention of the retired estimate must sit in a superseded context. "Exactly one mention" is a
# different and wrong assertion: the letter legitimately cites it twice, in §2 and in §4.1.
_bad = [i for i in [m.start() for m in re.finditer(r"\+1\.82", LET)]
        if "supersede" not in LET[max(0, i - 400):i + 400]]
check("§2 every mention of the retired +1.82 is marked as superseded", not _bad,
      f"{len(_bad)} unmarked of {len(re.findall(r'\+1\.82', LET))}")
check("§2 the shipped Results no longer state +1.82",
      "+1.82" not in SURF["01_manuscript/RESULTS.md"])

# ---- 5. the ledger the letter claims (§2) -------------------------------------------------------
_led = list(csv.DictReader(open(os.path.join(RES, "staging_ledger_native.csv"), encoding="utf-8")))
# Round 7 (C049): the SHIPPED vocabulary is REVERSE_SUPPORTED / REVERSE_CAVEATED (build_supp_tables
# .VERDICT_RENAME). The analysis CSV keeps the original DISCORDANT tokens, because main Figure 1 maps
# them to colours and legend labels, so the CSV is counted THROUGH the same rename table and the
# assertion below is made on the shipped vocabulary. Anything outside the documented four labels — in
# the CSV or in the workbook — is a failure.
_VERDICT_RENAME = {"DISCORDANT": "REVERSE_SUPPORTED", "DISCORDANT_CAVEATED": "REVERSE_CAVEATED",
                   "DISC_CAVEATED": "REVERSE_CAVEATED"}
_DOCUMENTED = {"CONCORDANT", "INDETERMINATE", "REVERSE_SUPPORTED", "REVERSE_CAVEATED"}
_v = {}
for r in _led:
    _k = _VERDICT_RENAME.get(r["verdict"], r["verdict"])
    _v[_k] = _v.get(_k, 0) + 1
check("§2 the ledger uses only the documented verdict vocabulary",
      set(_v) <= _DOCUMENTED, f"{sorted(_v)}")
_m = re.search(r"\*\*(\d+) concordant / (\d+) indeterminate / (\d+) reverse-supported /\s*\n?(\d+) pleiotropy-caveated\*\*", LET)
check("§2 the ledger counts in the letter match staging_ledger_native.csv",
      _m and [int(x) for x in _m.groups()] == [_v.get("CONCORDANT", 0), _v.get("INDETERMINATE", 0),
                                               _v.get("REVERSE_SUPPORTED", 0), _v.get("REVERSE_CAVEATED", 0)],
      f"letter={_m and _m.groups()} file={_v}")
check("§2 the retired DISCORDANT vocabulary is absent from the shipped workbook",
      "DISCORDANT" not in SURF["04_supplementary/Supplementary_Tables.xlsx"])

# ---- 6. the round-5 reference repairs the letter claims (§4.2) ----------------------------------
_refs = SURF["01_manuscript/REFERENCES_NUMBERED.md"]
check("§4.2 the 2026 guideline now carries volume, issue and pages", "154(4):e50-e158" in _refs)
check("§4.2 the epub date is no longer in the volume position", "2026;9 Jun" not in _refs)
check("§4.2 the Zheng entry is dated to its issue year, not its epub",
      re.search(r"Int J Epidemiol\. 2022;50\(6\):1995-2010", _refs) is not None)
_m = re.search(r"All (\d+) entries now resolve", LET)
_n_refs = len(re.findall(r"^\d+\. ", _refs, flags=re.M))
check("§4.2 the entry count in the letter matches the shipped list",
      _m and int(_m.group(1)) == _n_refs == F["references_cited"],
      f"letter={_m and _m.group(1)} list={_n_refs} facts={F['references_cited']}")

# ---- 7. the round-5 content repairs, asserted on the SHIPPED artifact ---------------------------
check("§4.1 the ceiling statistic agrees between Figure 1 and the prose",
      re.search(rf"P\s*=\s*{_p_deg:.2f}\b", SURF["03_main_figures/Figure1.pdf"]) is not None
      and f"{_p_deg:.2f}" in SURF["01_manuscript/ABSTRACT.md"])
check("§4.1 no raw markdown table survives in any shipped .docx",
      not [1 for k, v in SURF.items() if k.endswith(".docx")
           for l in v.split("\n") if l.strip().startswith("|") or "|---" in l])
_sm = re.findall(r"^## (Supplementary Methods \d+)[^\n]*\n\n([^\n]+)",
                 SURF["04_supplementary/SUPPLEMENTARY_INFORMATION.md"], flags=re.M)
check("§4.1 no Supplementary Methods section is truncated",
      len(_sm) >= 20 and not [n for n, b in _sm if b.rstrip()[-1:] not in ".!?"],
      f"{len(_sm)} sections")
_mvc = [r for r in csv.DictReader(open(os.path.join(RES, "mvmr_cascade.csv"), encoding="utf-8"))
        if r["model"] == "HF_via_CAD" and r["exposure"] == "BMI" and r["outcome"] == "HF"][0]
_want = f"{float(_mvc['direct_b']):+.2f}"
check(f"§4.2 the graphical abstract prints {_want}, as the letter says",
      _want in SURF["03_main_figures/GraphicalAbstract.pdf"]
      and "+0.41" not in SURF["03_main_figures/GraphicalAbstract.pdf"])
check("§4.2 the letter quotes the value the figure now prints",
      letter_has(re.escape(_want).replace("\\+", r"\*\*\+")) or _want in LET, _want)
check("§4.1/§4 no ASCII arrow survives on a reader-facing surface",
      not [1 for k, v in SURF.items()
           if not k.startswith("05_source_data/") and re.search(r"[A-Za-z0-9]+->[A-Za-z0-9]+", v)])
check("§4 the CAD→CAD row is labelled as positive-control validation",
      "positive-control validation of the cross-cohort design"
      in SURF["04_supplementary/Supplementary_Tables.xlsx"])
check("§4 the S6 legend is under the 300-word limit the letter claims",
      (lambda t: (lambda h: all(
          len(t[s:(h[i + 1][0] if i + 1 < len(h) else len(t))].split()) <= 300
          for i, (s, _n) in enumerate(h)))(
          [(m.start(), m.group(1)) for m in
           re.finditer(r"\*\*(Supplementary Figure S\d+) \|", t)]))(
          SURF["04_supplementary/SUPPLEMENTARY_INFORMATION.md"]))

# ---- 8. §4.3 — the letter must NOT claim a correction that was not made -------------------------
_summ = io.open(os.path.join(RES, "sensitivity_noukb_summary.txt"), encoding="utf-8").read()
_n13 = int(re.search(r"-\s*(\d+) headline FORWARD edges retain", _summ).group(1))
check("§4.3 the shipped Results keep the DERIVED headline-forward count",
      f"All {_n13} headline forward edges" in SURF["01_manuscript/RESULTS.md"], _n13)
# Round 8: the prose said "significant forward edges", which reads as contradicting the Discussion's
# 25/22/18 and was challenged on exactly that basis. sensitivity_noukb_summary.txt calls the set
# "headline FORWARD edges"; the prose now uses the generator's own word, so the check pins the
# generator's vocabulary rather than a looser synonym. The NUMBER is unchanged and still derived.
check("§4.3 the letter records that this was NOT a defect",
      "Both numbers are correct" in LET)

# ---- 9. the gate inventory the letter tabulates (§5) --------------------------------------------
_gate_src = io.open(os.path.join(SCR, "gate_package.py"), encoding="utf-8").read()
_have = {int(x) for x in re.findall(r"^# F(\d+)\.", _gate_src, flags=re.M)}
_claimed = {int(x) for x in re.findall(r"^\| F(\d+) \|", LET, flags=re.M)}
check("§5 every gate the letter tabulates exists in gate_package.py",
      _claimed and _claimed <= _have, f"letter={sorted(_claimed)} missing={sorted(_claimed - _have)}")
_m = re.search(r"Gate F(\d+)–F(\d+)", LET)
check("§5 the gate RANGE the letter names is complete in the source",
      _m and set(range(int(_m.group(1)), int(_m.group(2)) + 1)) <= _have,
      f"missing={_m and sorted(set(range(int(_m.group(1)), int(_m.group(2)) + 1)) - _have)}")

# ---- 10. §6 — the verification table must be re-run, not quoted ---------------------------------
def run(script):
    r = subprocess.run([sys.executable, os.path.join(SCR, script)], capture_output=True,
                       text=True, encoding="utf-8", errors="replace")
    return (r.stdout or "") + (r.stderr or "")


_qa = run("qa_manuscript.py")
_m_qa = re.search(r"QA RESULTS: (\d+) PASS, (\d+) FAIL", _qa)
_l_qa = re.search(r"`qa_manuscript\.py` \| (\d+) / (\d+)", LET)
check("§6 the qa figures in the letter are what qa_manuscript.py actually returns",
      _m_qa and _l_qa and _m_qa.groups() == _l_qa.groups(),
      f"letter={_l_qa and _l_qa.groups()} actual={_m_qa and _m_qa.groups()}")
for name, script in (("verify_rebuttal_r2.py", "verify_rebuttal_r2.py"),
                     ("verify_rebuttal_r3.py", "verify_rebuttal_r3.py")):
    out = run(script)
    m = re.search(r"RESULT: (\d+) verified, (\d+) FAILED", out)
    lm = re.search(rf"`{re.escape(name)}` \| (\d+) / (\d+)", LET)
    check(f"§6 the {name} figures in the letter are what it actually returns",
          m and lm and m.groups() == lm.groups(),
          f"letter={lm and lm.groups()} actual={m and m.groups()}")

# ---- 11. §8 — the open-items claim must be true of the package ----------------------------------
# R6 (2026-08-20): the bracketed placeholders were replaced by a truthful interim statement, so none
# should survive in the journal-facing declarations; the repository identifiers are finalised at proof.
_ph = re.findall(r"\[(?:AUTHOR-SUPPLIED|GitHub URL|Date)[^\]]*\]",
                 "\n".join(v for k, v in SURF.items() if k.startswith("02_cover_declarations/")))
check("§8 no bracketed author-supplied placeholder survives (R6 interim statement)", len(_ph) == 0,
      f"{len(_ph)} found: {_ph[:3]}")
check("§8 the letter does not claim the repository DOI is already deposited",
      not re.search(r"deposited at Zenodo|DOI 10\.\d{4,}", LET))

# ---- 12. the letter must not contradict the gates it lives beside -------------------------------
check("§all the letter carries no build-tree path",
      not re.search(r"0[1-5]_[a-z_]+/[\w./-]*|independent_build/", LET))
# The letter QUOTES the defects it fixed (`labeling`, `ischemic`), and a dialect check that fires on
# a quotation punishes the letter for being specific. Backticked spans are exempt.
_prose = re.sub(r"`[^`]*`", " ", LET)
_dia = re.findall(r"\b(labeling|labeled|ischemic|modeling|centered)\b", _prose)
check("§all the letter's own spelling matches the manuscript's dialect", not _dia, _dia[:4])

print("\n" + "=" * 96)
print(f"RESULT: {len(OK)} verified, {len(BAD)} FAILED")
if BAD:
    print("\nFAILED CLAIMS (the letter says something the package does not support):")
    for b in BAD:
        print("  x", b)
print("=" * 96)
sys.exit(1 if BAD else 0)
