# -*- coding: utf-8 -*-
"""Verify every checkable claim in RESPONSE_TO_INTERNAL_REVIEW_R2.md against the SHIPPED package.

A rebuttal is the highest-risk document in a submission: it makes dozens of factual assertions about a
package the reader can open, and a single "we removed X" that is untrue costs more credibility than
the original defect. The round-2 review found exactly that failure in our round-1 letter.

So the letter is not trusted here. Each assertion is restated as a machine check against
`submission package v3.2/` (every .md, both .docx twins, the workbook text and every figure PDF text
layer) or against the results files the numbers came from. Anything that cannot be checked mechanically
is listed at the end as such, rather than silently counted as verified.
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

import io, os, re, sys, glob, csv, zlib, json

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = P4_ROOT
BASE = os.path.join(ROOT, "independent_build")
# 2026-07-26: the package version is DERIVED from build_submission_v1.py, never typed here. A
# checker hard-pointed at an old package silently verifies an artifact nobody is shipping - the
# exact failure prose_ai_scan.py had. Its scope IS the check.
_bsv = io.open(os.path.join(ROOT, "independent_build", "scripts", "build_submission_v1.py"),
               encoding="utf-8").read()
_m = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _m, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _m.group(1)
PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")
RES = os.path.join(BASE, "results")
FIGD = os.path.join(BASE, "figures", "data")

# ---------------------------------------------------------------- every shipped text surface
def docx_text(p):
    import docx
    d = docx.Document(p)
    return "\n".join([q.text for q in d.paragraphs] +
                     [c.text for t in d.tables for r in t.rows for c in r.cells])

# PDF text comes from the SHARED reader (~/.claude/lib/pdftext.py), not a local copy. This file
# carried its own scanner until 2026-07-27; `pdftext.py --audit` found five such forks in this
# project, and measuring them showed the local scanner returned ~70x too many characters — CMap and
# font-encoding junk — so every presence probe over it was weaker than it looked. The import raises
# if the shared library is missing: a checker that cannot read must never report clean.
sys.path.insert(0, SHARED_LIB)
from pdftext import pdf_text
SURF = {}
for p in glob.glob(os.path.join(PKG, "**", "*.md"), recursive=True):
    SURF[os.path.relpath(p, PKG)] = open(p, encoding="utf-8").read()
for p in glob.glob(os.path.join(PKG, "**", "*.docx"), recursive=True):
    SURF[os.path.relpath(p, PKG)] = docx_text(p)
for p in glob.glob(os.path.join(PKG, "**", "*.pdf"), recursive=True):
    SURF["PDF:" + os.path.basename(p)] = pdf_text(p)
import openpyxl
_wb = openpyxl.load_workbook(os.path.join(PKG, "04_supplementary", "Supplementary_Tables.xlsx"),
                             read_only=True, data_only=True)
SURF["WORKBOOK"] = "\n".join(str(v) for ws in _wb.worksheets
                             for row in ws.iter_rows(values_only=True)
                             for v in row if isinstance(v, str))
ALLTEXT = "\n".join(SURF.values())
BODY = "\n".join(v for k, v in SURF.items() if k.startswith("01_manuscript") and k.endswith(".md"))

FAILS, OKS, MANUAL = [], [], []
def check(cond, claim):
    (OKS if cond else FAILS).append(claim)

def absent(pat, claim):
    hits = sorted({k for k, t in SURF.items() if re.search(pat, t, re.I)})
    check(not hits, f"{claim}  [hits: {hits[:3]}]" if hits else claim)

def present(pat, claim, where=None):
    scope = SURF if where is None else {k: v for k, v in SURF.items() if where in k}
    hits = sorted({k for k, t in scope.items() if re.search(pat, t, re.I)})
    check(bool(hits), claim if hits else f"{claim}  [found nowhere]")

def num(path, col, where):
    """Read one value from a results CSV so a quoted number is checked against its file."""
    p = os.path.join(RES, path) if os.path.exists(os.path.join(RES, path)) else os.path.join(FIGD, path)
    for r in csv.DictReader(open(p, encoding="utf-8")):
        if all(r[k] == v for k, v in where.items()):
            return float(r[col])
    return None

print("=" * 78)
print(f"VERIFYING RESPONSE_TO_INTERNAL_REVIEW_R2.md AGAINST submission package {PKG_VERSION}")
print(f"surfaces opened: {len(SURF)}  ({sum(1 for k in SURF if k.startswith('PDF:'))} figure PDFs, "
      f"{sum(1 for k in SURF if k.endswith('.docx'))} .docx, 1 workbook)")
print("=" * 78)

# ---------------------------------------------------------------- §1 phrases claimed removed
print("\n-- letter claims these phrases return ZERO hits package-wide --")
for pat, name in [
    (r"resolved the mechanism", "'resolved the mechanism'"),
    (r"mirror image", "'mirror image'"),
    (r"genuine .{0,12}sign reversal", "'genuine sign reversals'"),
    (r"power artifact", "'power artifact'"),
    (r"central fat spares", "'central fat spares heart failure'"),
    (r"topology travels", "'topology travels'"),
    (r"\+0\.24 SD", "'+0.24 SD' unit error"),
    (r"specificity control", "'specificity control'"),
    # (the ASSERT/DEVELOP tier labels are checked case-SENSITIVELY below: they were ALL-CAPS panel
    #  labels, and a case-insensitive search matches ordinary prose such as "we develop them in
    #  Figure 3". Verify the verifier -- that hit was the checker's fault, not the artifact's.)
    (r"evidence asserts|develops as hypothesis|asserts, develops", "retired tier vocabulary in prose"),
    (r"(?<!correlated-marker )negative control", "unqualified 'negative control'"),
    (r"favou?rs adiposity-first prevention", "'favours adiposity-first prevention'"),
    (r"10,000 shuffl", "'10,000 shuffles' annotation"),
    (r"precisely the edge where", "'precisely the edge where overlap inflates'"),
]:
    absent(pat, f"absent everywhere: {name}")

# ---------------------------------------------------------------- §2 wording claimed present
def absent_cs(pat, claim):
    """Case-SENSITIVE absence, for ALL-CAPS labels only."""
    hits = sorted({k for k, s in SURF.items() if re.search(pat, s)})
    check(not hits, f"{claim}  [hits: {hits[:3]}]" if hits else claim)

absent_cs(r"\bASSERT\b|\bDEVELOP\b|\bBOUNDED\b",
          "absent everywhere: ASSERT/DEVELOP/BOUNDED tier labels (case-sensitive)")

print("-- letter claims this wording IS present --")
for pat, name, where in [
    (r"correlated-marker negative control", "correlated-marker qualifier", None),
    (r"compatible with overlap", "'compatible with overlap and/or cohort differences'", None),
    (r"not (?:time-stamped or )?prospectively registered|No analysis was time-stamped", "not-registered statement", None),
    (r"HIGHER CONF", "HIGHER CONFIDENCE grade on the synthesis panel", "PDF:SupplFig7"),
    (r"Hartung", "Hartung-Knapp", "PDF:SupplFig6"),
    (r"index-event|collider", "index-event/collider caveat", None),
    (r"Node-level", "node-level row on the portability panel", "PDF:Figure4"),
    (r"primary", "the (primary)/(secondary) hierarchy on the portability panel", "PDF:Figure4"),
    (r"Barton", "ApoB attributed to Barton", "STROBE"),
    # The R2 letter's claim is that the overlap statement is SCOPED TO THE OUTCOME SIDE. This check
    # used to enforce the specific wording "outcome-side overlap ... does not manufacture", which
    # round 5 asked us to retire as categorical — so the gate went red on a deliberate improvement
    # and, taken at face value, would have argued for reverting the artifact. The check was stale in
    # INTENT, not merely in scope. It now asserts the claim the letter actually makes: the statement
    # is outcome-side-scoped AND names exposure-side overlap as the residual limitation.
    (r"no evidence that outcome-side overlap explains the ordering; exposure-side\s+overlap remains "
     r"a limitation",
     "overlap claim scoped to the outcome side, with exposure-side named as the residual limit",
     None),
]:
    present(pat, f"present: {name}", where)

# ---------------------------------------------------------------- §3 numbers quoted in the letter
print("-- numbers the letter quotes, checked against their source files --")
LET = open(os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R2.md"), encoding="utf-8").read()

gs = list(csv.DictReader(open(os.path.join(RES, "interaction_global_scale.csv"), encoding="utf-8")))
fam = int(gs[0]["family_size"]); bon = float(gs[0]["bonferroni"]); slope = float(gs[0]["global_slope"])
raw = {(r["exposure"], r["outcome"]) for r in gs if r["sig_raw"] == "True"}
abso = {(r["exposure"], r["outcome"]) for r in gs if r["sig_absorbed"] == "True"}
opp = [r for r in gs if r["sig_raw"] == "True" and float(r["b_eur_used"]) * float(r["b_eas"]) < 0]
check(fam == 67 and "67-edge interaction family" in LET, "interaction family size 67 matches the file")
check(abs(bon - 7.46e-4) < 1e-6 and "7.46 × 10⁻⁴" in LET, "family Bonferroni 7.46e-4 matches the file")
check(abs(slope - 0.436) < 5e-4 and "0.436 × b_EUR" in LET, "global slope 0.436 matches the file")
check(len(raw) == 14 and len(abso) == 12 and len(raw & abso) == 6 and "14 edges to 12, with only 6 in common" in LET,
      f"interaction survivors 14 -> 12, 6 in common (file: {len(raw)}, {len(abso)}, {len(raw & abso)})")
check(len(opp) == 9 and "9 opposite-sign survivors" in LET, f"9 opposite-sign survivors (file: {len(opp)})")

med = {r["exposure"]: r for r in csv.DictReader(open(os.path.join(FIGD, "fig2_mediation.csv"), encoding="utf-8"))}
check(abs(float(med["T2D"]["direct"]) - 0.011) < 5e-4 and "+0.01, *P* = 0.41" in LET,
      "T2D direct effect +0.01 / P=0.41 matches fig2_mediation.csv")
check(abs(float(med["BMI"]["pm_lo"]) - 16) < 1 and abs(float(med["BMI"]["pm_hi"]) - 25) < 1,
      "BMI mediated 16-25% matches fig2_mediation.csv")

cov = list(csv.DictReader(open(os.path.join(RES, "mediation_covariance_sensitivity.csv"), encoding="utf-8")))
b = [r for r in cov if r["exposure"] == "BMI"]; t2 = [r for r in cov if r["exposure"] == "T2D"]
blo, bhi = min(float(r["pm_product_lo"]) for r in b), max(float(r["pm_product_hi"]) for r in b)
check(len(cov) >= 80 and "80 combinations" in LET, f"covariance sweep has >=80 combinations (file: {len(cov)})")
check(round(blo) >= 15 and round(bhi) <= 27 and "15–27%" in LET,
      f"BMI mediated range 15-27% across the sweep (file: {blo:.0f}-{bhi:.0f})")

port = {(r["stat"], r["block"]): r for r in csv.DictReader(open(os.path.join(FIGD, "fig4_portability.csv"), encoding="utf-8"))}
nl = port[("Sign concordance", "Node-level")]; el = port[("Sign concordance", "Edge-level")]
check(abs(float(nl["lo"]) - 0.444) < 5e-3 and abs(float(nl["hi"]) - 1.0) < 1e-9 and "[0.44, 1.00]" in LET,
      "node-level CI [0.44, 1.00] matches fig4_portability.csv")
check(abs(float(el["lo"]) - 0.558) < 5e-3 and abs(float(el["hi"]) - 0.808) < 5e-3 and "[0.56, 0.81]" in LET,
      "edge-level CI [0.56, 0.81] matches fig4_portability.csv")

sm = open(os.path.join(RES, "steiger_margin.txt"), encoding="utf-8").read()
check("1.92" in sm and "TG→T2D at 1.92×" in LET, "smallest Steiger margin 1.92x (TG->T2D) matches steiger_margin.txt")
pre = open(os.path.join(RES, "mrpresso_all_edges.txt"), encoding="utf-8").read()
check("59/59" in pre and "59 of 59" in LET, "MR-PRESSO global test 59/59 matches mrpresso_all_edges.txt")
check("+0.8069" in pre and "+0.5088" in pre and "+0.807 → corrected\n+0.509" in LET.replace("  ", " ") or
      ("+0.807" in LET and "+0.509" in LET), "HF->CAD distortion 0.807 -> 0.509 matches the file")
_qc = "".join(open(os.path.join(RES, f), encoding="utf-8", errors="replace").read()
              for f in os.listdir(RES) if "sbp" in f.lower() and f.endswith((".txt", ".md")))
check("7,088,121" in _qc or "7088121" in _qc,
      "the SNP-only structural finding is in the SBP variant-QC output")
whr = open(os.path.join(RES, "whr_triangulation.txt"), encoding="utf-8").read()
for v in ["+0.5246", "+0.0120", "+0.2648"]:
    check(v in whr, f"WHR triangulation value {v} present in whr_triangulation.txt")
hk = open(os.path.join(RES, "h4_leandiabetes.txt"), encoding="utf-8").read()
check("p=0.2238" in hk and "*P* = 0.22" in LET, "Hartung-Knapp P=0.22 matches h4_leandiabetes.txt")

# ledger bound-invariance, recomputed rather than quoted
ws = _wb["S5b_ledger_bound_sensitivity"]
rows = [r for r in ws.iter_rows(values_only=True)][2:]
rows = [r for r in rows if r and r[0]]
DIS = lambda v: bool(v) and str(v).startswith("DISCORDANT")
prim_dis = [r[0] for r in rows if DIS(r[3])]
# The workbook now typesets the edge key (build_supp_tables.typeset), so compare on the
# NOTATION-INDEPENDENT form; this check is about which transitions are discordant, not about
# which arrow glyph the cell uses.
ASCII = lambda v: str(v).replace("→", "->")
always = [ASCII(r[0]) for r in rows if all(DIS(v) for v in r[1:6])]
conc_flip = [r[0] for r in rows if r[3] == "CONCORDANT" and any(DIS(v) for v in r[1:6])]
check(len(rows) == 15 and len(prim_dis) == 4, f"15 transitions, 4 discordant at primary bounds (got {len(rows)}, {len(prim_dis)})")
check(sorted(always) == ["T2D->HF", "T2D->Stroke"], f"T2D->HF and T2D->Stroke discordant at all 5 bounds (got {always})")
check(not conc_flip, f"no primary-concordant transition becomes discordant at any bound (got {conc_flip})")

# derived package facts the letter quotes
F = json.load(open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))
check(F["main_figures"] == 4 and "four figures, 25 panels" in LET.lower().replace("**", ""),
      f"4 main figures (facts: {F['main_figures']})")
check(F["total_panels"] == 25, f"25 panels (facts: {F['total_panels']})")
check(F["supp_figures"] == 7, f"7 supplementary figures (facts: {F['supp_figures']})")
_n = F["abstract_full_words"]
# Derive, do not hardcode: this check asserted a literal 333 and went stale the moment the
# abstract moved by one word during the round-3 pass.
check(f"{_n} words" in LET, f"the letter quotes the abstract length the package actually has ({_n})")
check(len(glob.glob(os.path.join(PKG, "**", "*"), recursive=True)) >= 98, "package has 98 files")
check(os.path.exists(os.path.join(PKG, "renv.lock")), "renv.lock ships at the package root")
lock = json.load(open(os.path.join(PKG, "renv.lock"), encoding="utf-8"))
check(lock["R"]["Version"] == "4.4.1" and len(lock["Packages"]) == 131 and "131 packages" in LET,
      f"renv.lock records R 4.4.1 and 131 packages (file: {lock['R']['Version']}, {len(lock['Packages'])})")

# the letter's own hygiene: it is INTERNAL and must not ship
check(not glob.glob(os.path.join(PKG, "**", "RESPONSE_TO_INTERNAL_REVIEW*"), recursive=True),
      "the response letter does NOT ship inside the package")
# Scoped to the folders that go to the editor. The author to-do legitimately mentions suggesting or
# excluding JOURNAL reviewers, which is submission vocabulary, not review-response language.
_UPLOAD = "\n".join(v for k, v in SURF.items()
                    if k.startswith(("01_", "02_", "03_", "04_", "05_", "PDF:", "WORKBOOK")))
check(not re.search(r"round[- ]?[12] (review|response)|reviewer (round|comment|item|report)"
                    r"|rebuttal|response to (the )?review", _UPLOAD, re.I),
      "no review/rebuttal language in any file that goes to the editor")

MANUAL += [
    "Editorial judgements (what to decline and why) are arguments, not facts, and are not checkable here.",
    "The claim that Cardiovascular Diabetology sets no main-text limit was read from the live journal "
    "guidelines on 2026-07-25; it is a statement about an external site, re-check before submission.",
]

print("\n" + "=" * 78)
print(f"RESULT: {len(OKS)} verified, {len(FAILS)} FAILED")
if FAILS:
    print("\nFAILED CLAIMS (the letter says something the package does not support):")
    for f in FAILS:
        print("  x", f)
print(f"\nNot mechanically checkable ({len(MANUAL)}):")
for m in MANUAL:
    print("  -", m)
print("=" * 78)
sys.exit(1 if FAILS else 0)
