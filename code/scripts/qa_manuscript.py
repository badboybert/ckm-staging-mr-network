# -*- coding: utf-8 -*-
"""Pre-submission QA for CKM Paper 4: file existence, figure/table cross-reference callouts,
citation<->reference correspondence, MR causal-language / calibration-invariant compliance, and
key-number presence across sections. Deterministic; prints PASS/FAIL per check + a final tally."""
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

import os, re, glob

BASE = P4_BASE
FIG  = os.path.join(BASE, "figures")
RES  = os.path.join(BASE, "results")
MAN  = os.path.join(BASE, "manuscript")

def rd(p):
    return open(p, encoding="utf-8").read() if os.path.exists(p) else ""

def shipped(txt):
    """Strip HTML provenance comments before any check.

    build_manuscript_docx.py strips these comments, so they never reach the reader. Scanning them
    made the gate fail on its own changelog -- e.g. the RESULTS header comment records that this
    draft 'supersedes the 111-edge draft', which is a true statement ABOUT the file, not a stale
    claim IN the manuscript. A gate must assert on what ships."""
    return re.sub(r"<!--.*?-->", "", txt, flags=re.S)

SEC = {n: shipped(rd(os.path.join(FIG, f"{n}.md"))) for n in
       ["INTRODUCTION","METHODS","RESULTS","DISCUSSION","CONCLUSIONS","ABBREVIATIONS",
        "ABSTRACT","FIGURE_LEGENDS","SUPPL_FIGURES"]}
BODY = (SEC["INTRODUCTION"] + SEC["METHODS"] + SEC["RESULTS"] + SEC["DISCUSSION"]
        + SEC["CONCLUSIONS"])

fails, warns, oks = [], [], []
def check(cond, msg):
    (oks if cond else fails).append(msg)
def warn(cond, msg):
    if not cond: warns.append(msg)

# ---------- 1. EXISTENCE ----------
# The figure ranges are DERIVED from the legend files rather than written as literals. Hardcoded
# ranges (1..7 main, 1..6 supplementary) had to be edited by hand every time the main/supplementary
# split moved, and on 2026-07-25 (Figures 5 and 6 -> Supplementary S6 and S7) they failed for the
# wrong reason: they asserted the CONTINUED EXISTENCE of files the split had deliberately removed.
_main_ids = sorted(int(x) for x in re.findall(r"^\*\*Figure (\d+) \|",
                                              SEC["FIGURE_LEGENDS"], flags=re.M))
_supp_ids = sorted(int(x) for x in re.findall(r"^\*\*Supplementary Figure S(\d+) \|",
                                              SEC["SUPPL_FIGURES"], flags=re.M))
check(_main_ids == list(range(1, len(_main_ids) + 1)),
      f"main figure legends are numbered contiguously from 1 (found {_main_ids})")
check(_supp_ids == list(range(1, len(_supp_ids) + 1)),
      f"supplementary figure legends are numbered contiguously from 1 (found {_supp_ids})")
must_exist = (
    [os.path.join(FIG, f"Figure{i}.{e}") for i in _main_ids for e in ("png", "pdf")] +
    [os.path.join(FIG, f"SupplFig{i}.{e}") for i in _supp_ids for e in ("png", "pdf")] +
    [os.path.join(FIG, f"{n}.md") for n in ["INTRODUCTION","METHODS","RESULTS","DISCUSSION",
        "CONCLUSIONS","ABBREVIATIONS",
        "ABSTRACT","FIGURE_LEGENDS","SUPPL_FIGURES","REFERENCES_MASTER","REFERENCES_NUMBERED","FRONT_MATTER"]] +
    [os.path.join(MAN, "CKM_Paper4_manuscript.docx"), os.path.join(MAN, "Supplementary_Tables.xlsx")]
)
# A retired figure file must not linger in the tree either: it would be copied by a future build.
for _stray in sorted(glob.glob(os.path.join(FIG, "Figure*.pdf"))):
    _n = int(re.search(r"Figure(\d+)\.pdf$", os.path.basename(_stray)).group(1))
    check(_n in _main_ids, f"NO STRAY main-figure file: Figure{_n}.pdf has no legend")
for p in must_exist:
    check(os.path.exists(p), f"EXISTS: {os.path.relpath(p, BASE)}")

# ---------- 2. MAIN-FIGURE CALLOUTS vs legends ----------
# legends define "Figure N | ..." and panels "(a)".."(e)"
leg = SEC["FIGURE_LEGENDS"]
panels_def = {}
for m in re.finditer(r"\*\*Figure (\d) \|", leg):
    n = m.group(1)
    seg = leg[m.end(): leg.find("**Figure", m.end()) if leg.find("**Figure", m.end())>0 else len(leg)]
    panels_def[n] = set(re.findall(r"\(\*{0,2}(\w)\*{0,2}\)", seg))   # legends use bold panels (**a**)
# Do not hardcode the figure count: assert the legends and the rendered files agree. A literal here
# would have to be edited every time the main/supplementary split moves (it did, 6 -> 4 on
# 2026-07-25), and a stale literal fails for the wrong reason or passes for none.
_fig_pdfs = sorted(int(m.group(1)) for m in
                   (re.match(r"Figure(\d+)\.pdf$", os.path.basename(p))
                    for p in glob.glob(os.path.join(FIG, "Figure*.pdf"))) if m)
check(sorted(int(k) for k in panels_def) == _fig_pdfs,
      f"FIGURE_LEGENDS defines the same figures that are rendered "
      f"(legends {sorted(int(k) for k in panels_def)} vs files {_fig_pdfs})")
# callouts in Results/Discussion like "Figure 1e", "Figure 2a", "Figure 4b,c", "Figure 5a–c"
callout_txt = SEC["RESULTS"] + SEC["DISCUSSION"]
bad_callouts = []
for m in re.finditer(r"Figure (\d)([a-z](?:[,–\-][a-z])*)?", callout_txt):
    fn, pnls = m.group(1), m.group(2)
    if fn not in panels_def:
        bad_callouts.append(f"Figure {fn} (no legend)"); continue
    if pnls:
        for pl in re.findall(r"[a-z]", pnls):
            if pl not in panels_def[fn]:
                bad_callouts.append(f"Figure {fn}{pl} (panel not in legend {sorted(panels_def[fn])})")
check(not bad_callouts, f"main-figure panel callouts all defined in legends" + (f" -- BAD: {bad_callouts}" if bad_callouts else ""))

# ---------- 3. SUPPLEMENTARY callouts exist ----------
# Supp figures S1-S6 (SUPPL_FIGURES.md), Supp tables S1-S12 (workbook)
supp_fig_callouts = set(re.findall(r"Supplementary Fig(?:ure|\.)?\s*S?(\d+)", BODY, re.I))
supp_tab_callouts = set(re.findall(r"Supplementary Table\s*S?(\d+)", BODY, re.I))
sf_defined = set(re.findall(r"Supplementary Figure S(\d)", SEC["SUPPL_FIGURES"]))
for n in sorted(supp_fig_callouts):
    check(n in sf_defined, f"Supplementary Fig S{n} callout -> legend exists")
# workbook sheet names
try:
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(MAN, "Supplementary_Tables.xlsx"), read_only=True)
    sheet_tabs = set(re.findall(r"S(\d+)", " ".join(wb.sheetnames)))
except Exception as e:
    sheet_tabs = set(); warns.append(f"could not read workbook: {e}")
for n in sorted(supp_tab_callouts):
    check(n in sheet_tabs, f"Supplementary Table S{n} callout -> workbook sheet exists")

# ---------- 4. CITATION <-> REFERENCE correspondence ----------
numbered = rd(os.path.join(FIG, "REFERENCES_NUMBERED.md"))
m = re.search(r"(\d+) references cited", numbered)
check(bool(m), "REFERENCES_NUMBERED.md reports a cited count")
if m:
    ncited = int(m.group(1))
    check("0 uncited" not in numbered or "uncited" not in numbered.lower() or True, "")  # placeholder
# re-run the numbered builder's core logic by importing counts from its output tail is unreliable;
# instead re-execute it and capture stdout
import subprocess, sys
r = subprocess.run([sys.executable, os.path.join(BASE,"scripts","build_numbered_refs.py")],
                   capture_output=True, text=True)
out = r.stdout
mm = re.search(r"(\d+) cited, (\d+) uncited, (\d+) unmatched", out)
if mm:
    c,u,x = map(int, mm.groups())
    check(u==0, f"references: {u} uncited (want 0)")
    check(x==0, f"references: {x} unmatched (want 0)")
    check(c>=45, f"references cited = {c} (>=45)")
else:
    fails.append("could not parse build_numbered_refs output")

# ---------- 5. MR causal-language / forbidden phrases ----------
FORBIDDEN = [r"proves\s+caus", r"demonstrat\w*\s+caus", r"establish\w*\s+that\s+\w+\s+caus",
             r"\bproven\s+to\s+cause\b", r"definitively\s+caus"]
for pat in FORBIDDEN:
    hits = [n for n in ["INTRODUCTION","METHODS","RESULTS","DISCUSSION","ABSTRACT"] if re.search(pat, SEC[n], re.I)]
    check(not hits, f"no forbidden causal phrase /{pat}/" + (f" -- in {hits}" if hits else ""))
# banned staging mis-rounding: the PRIMARY 0.926 concordance must pair with 1.2e-3 (the 132-edge
# value), never with the superseded 1.8e-3 (111-edge) or a bare 1e-3.
# (UKB-free staging is legitimately 1.1e-3, paired with concordance 1.000 -- so check proximity to 0.926.)
for n in ["RESULTS","DISCUSSION","FIGURE_LEGENDS","ABSTRACT"]:
    for mm in re.finditer(r"0\.926", SEC[n]):
        window = SEC[n][mm.start(): mm.start()+90]
        check("1.8 × 10⁻³" not in window, f"{n}: 0.926 not mis-paired with the superseded '1.8 × 10⁻³'")

# superseded values must not reappear anywhere in the main text (the 111-edge network and the
# retired Kendall-tau / binomial portability statistics).
SUPERSEDED = {
    # NB: a bare \b111\b is WRONG here -- the Results legitimately says the staging result "survived
    # completing the network from 111 to 132 edges", a deliberate historical reference. Only stale
    # USES of 111 as the current network size are banned.
    "111 as current network size": r"111[- ]edge|111 directed|network of 111|all 111\b",
    "0.05/111 Bonferroni":   r"0\.05/111",
    "111-edge staging P":    r"1\.8 × 10⁻³",
    "retired Kendall tau":   r"Kendall|τ = 0\.49",
    "retired binomial 44/65": r"44 of 65|44/65",
    "abandoned exclusive T2D": r"only via CAD|reaches heart failure only",
    "pre-registered (invariant #10)": r"[Pp]re-registered",
}
for label, pat in SUPERSEDED.items():
    hits = [n for n in ["INTRODUCTION","METHODS","RESULTS","DISCUSSION","ABSTRACT","FIGURE_LEGENDS"]
            if re.search(pat, SEC[n])]
    check(not hits, f"superseded '{label}' absent" + (f" -- FOUND in {hits}" if hits else ""))

# ---------- 6. KEY-NUMBER presence across sections ----------
KEY = {
  "CAD->HF 0.285":        (["RESULTS","DISCUSSION","ABSTRACT"], "0.285"),
  "CAD->HF P 1e-119":     (["RESULTS","DISCUSSION","ABSTRACT"], "10⁻¹¹⁹"),
  "staging conc 0.926":   (["RESULTS","DISCUSSION","FIGURE_LEGENDS"], "0.926"),
  "staging P 1.5e-3 (exact)":(["RESULTS","DISCUSSION","FIGURE_LEGENDS","ABSTRACT","METHODS"], "1.5 × 10⁻³"),
  # Scope narrowed to RESULTS on 2026-07-25: the Cardiovascular Diabetology 350-word abstract limit
  # cannot carry the reverse-edge Steiger P. What mattered about it in the abstract was the CALIBRATION
  # -- that CAD->HF must not read as a resolved direction -- so the number check is replaced there by
  # a claim check on the hedge itself, asserted immediately below. That is a stronger test, not a
  # weaker one: the number could be present while the hedge was dropped.
  "HF->CAD Steiger 2e-78":(["RESULTS"], "10⁻⁷⁸"),
  "SBP->CKD EUR 1e-17":   (["RESULTS","DISCUSSION","ABSTRACT"], "10⁻¹⁷"),
  "SBP->CKD EAS 0.66":    (["RESULTS","DISCUSSION","ABSTRACT"], "0.66"),
  "BMI->HF 1e-31":        (["RESULTS","DISCUSSION","ABSTRACT"], "10⁻³¹"),
  "T2D->HF direct 0.41":  (["RESULTS","DISCUSSION","ABSTRACT"], "0.41"),
  "132 edges":            (["RESULTS","DISCUSSION","METHODS","ABSTRACT","INTRODUCTION"], "132"),
  # UKB-free staging was re-run under EXACT enumeration on 2026-07-24. 3 × 10⁻⁴ was a 10,000-draw
  # Monte-Carlo estimate and is not attainable on the exact support (min non-zero 1/1980 = 5.1 × 10⁻⁴).
  "UKBfree staging exact P": (["RESULTS"], "7.6 × 10⁻⁴"),
  # --- the six recalibrations of the 132-edge rewrite must each be visible in the main text ---
  # 0.31 -> 0.27 on 2026-07-26: the shipped value came from a chain of 162 ATTEMPTS at 8.9%
  # acceptance, i.e. ~0.27 successful swaps per edge. The mixed chain (50 successful swaps/edge)
  # gives 0.27, flat from 10x to 100x across three seeds. Conclusion unchanged; the number moved.
  "degree-null 0.27":     (["RESULTS","DISCUSSION","ABSTRACT"], "0.27"),
  # Mediation is now reported as covariance-dependent ranges (round-2 #1), not a firm 68-79%.
  "mediation covariance range": (["RESULTS","FIGURE_LEGENDS"], "44–116"),
  "mediation CI 16-25%":  (["RESULTS","FIGURE_LEGENDS"], "16–25"),
  "CAC staging 0.941":    (["RESULTS","FIGURE_LEGENDS"], "0.941"),
  "interaction ns 0.43":  (["RESULTS","DISCUSSION","ABSTRACT"], "0.43"),
  # Was "1.1 × 10⁻³" over ["RESULTS","METHODS"] — the 111-edge value, which this gate went on
  # asserting after the network was completed. The UKB-free staging is now rebuilt on 132 edges.
  "UKBfree staging MAIN":  (["RESULTS"], "0.960"),
  "UKBfree staging FinnGen":(["RESULTS"], "0.913"),
}
for label,(files,val) in KEY.items():
    missing = [n for n in files if val not in SEC[n]]
    check(not missing, f"key-number '{label}' ({val}) present in {files}" + (f" -- MISSING in {missing}" if missing else ""))

# ---------- 6b. CALIBRATION CLAIMS, not just numbers ----------
# A number can survive while the hedge that makes it honest is trimmed away. These assert the hedge.
CLAIMS = {
    "abstract hedges the CAD->HF direction as unresolved":
        ("ABSTRACT", r"reverse direction unresolved|reverse .{0,40}(unresolved|not supported|directionally suspect)"),
    "abstract states the degree-null does not corroborate beyond trait roles":
        ("ABSTRACT", r"degree-preserving .{0,40}not exceeded"),
    "abstract keeps portability suggestive, not established":
        ("ABSTRACT", r"suggestive"),
    "abstract keeps the T2D mediated proportion imprecise":
        # The check asserts the HEDGE, not one phrasing of it. Round 4 reworded this to
        # "imprecisely quantified"; the old literal would have failed a sentence that still carries
        # the caveat, which is a gate testing its own memory rather than the claim.
        ("ABSTRACT", r"imprecise and covariance-dependent|imprecisely (quantified|estimated)"),
    "results reports the node-level interval as the primary one":
        ("RESULTS", r"dependence-robust interval as the primary"),
}
for label, (sec, pat) in CLAIMS.items():
    check(bool(re.search(pat, SEC[sec], re.I)), f"calibration claim -- {label}")

# ---------- report ----------
print(f"\n===== QA RESULTS: {len(oks)} PASS, {len(fails)} FAIL, {len(warns)} WARN =====\n")
if fails:
    print("FAILURES:")
    for f in fails: print("  ✗", f)
if warns:
    print("\nWARNINGS:")
    for w in warns: print("  !", w)
print(f"\n(passed {len(oks)} checks)")
