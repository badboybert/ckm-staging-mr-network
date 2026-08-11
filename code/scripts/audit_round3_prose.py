# -*- coding: utf-8 -*-
"""Row-by-row audit of the round-3 review's SECTION 10 prose revisions (and 6.3 diction).

The review's section 10 is a table of ~65 "original -> recommended" prose rows. Earlier sessions
closed the review's *analysis* items and its 6.3 diction list, but section 10 was never enumerated
row by row, so "what did we do with the prose suggestions" had no mechanical answer.

Each row is restated as two searches over every shipped surface (.md, both .docx twins, the workbook
cell text, and the text layer of all 12 figure PDFs):
  ORIGINAL  - a distinctive fragment of the wording the review quotes as shipped.
  ADOPTED   - a distinctive fragment of the substance the review recommends.

Verdict:
  DONE        original absent, recommended substance present
  DONE-DIFF   original absent, recommended fragment absent (fixed with different wording -> read it)
  NOT DONE    original still present
  DECLINED    recorded as a deliberate decline (still reports whether the original is present)

A failing row is a HYPOTHESIS. Read the artifact before acting on it.
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

import io, os, re, sys, glob, zlib

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = P4_ROOT
# 2026-07-26: the package version is DERIVED from build_submission_v1.py, never typed here. A
# checker hard-pointed at an old package silently verifies an artifact nobody is shipping - the
# exact failure prose_ai_scan.py had. Its scope IS the check.
_bsv = io.open(os.path.join(ROOT, "independent_build", "scripts", "build_submission_v1.py"),
               encoding="utf-8").read()
_m = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _m, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _m.group(1)
PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")


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
    SURF[os.path.relpath(p, PKG).replace("\\", "/")] = open(p, encoding="utf-8").read()
for p in glob.glob(os.path.join(PKG, "**", "*.docx"), recursive=True):
    SURF[os.path.relpath(p, PKG).replace("\\", "/")] = docx_text(p)
for p in glob.glob(os.path.join(PKG, "**", "*.pdf"), recursive=True):
    SURF["PDF:" + os.path.basename(p)] = pdf_text(p)
import openpyxl
_wb = openpyxl.load_workbook(os.path.join(PKG, "04_supplementary", "Supplementary_Tables.xlsx"),
                             read_only=True, data_only=True)
SURF["WORKBOOK"] = "\n".join(str(v) for ws in _wb.worksheets
                             for row in ws.iter_rows(values_only=True)
                             for v in row if isinstance(v, str))

assert len(SURF) >= 35, f"only {len(SURF)} surfaces loaded"

# ---------------------------------------------------------------- rows
# (id, section, original_regex|None, adopted_regex|None, status_override|None, note)
ROWS = [
    # ---- 10.1 Title, abstract, plain-language summary
    ("10.1-1", "Title", None, None, "DECLINED",
     "Title change to R3 wording; author chose the current title 2026-07-25."),
    ("10.1-2", "Abstract Background",
     r"is genetically untested", r"sequential mediation|common upstream driver", None, ""),
    ("10.1-3", "Abstract Methods",
     r"assessed its portability", r"132 (directed|non-self)", None, ""),
    ("10.1-4", "Abstract Staging",
     r"corroborating the order", r"degree-preserving", None, ""),
    ("10.1-5", "Abstract BMI",
     r"retained direct effects", r"CAD-independent", None, ""),
    ("10.1-6", "Abstract T2D",
     r"predominantly, not exclusively", r"no CAD-independent|not detected", None, ""),
    ("10.1-7", "Abstract Conclusion",
     r"direct driver of heart failure", r"beyond (coronary|CAD)", None, ""),
    ("10.1-8", "Plain-language summary",
     r"more than high blood sugar", None, None,
     "CVD has no plain-language-summary item; check whether the surface exists at all."),
    # ---- 10.2 Introduction
    ("10.2-1", "Intro p1", r"names what clinicians have long managed", None, None, ""),
    ("10.2-2", "Intro p1", r"imposed an order on it", r"clinical stages|staged clinical", None, ""),
    ("10.2-3", "Intro p1", r"arrows in fact point the way", None, None, ""),
    ("10.2-4", "Intro p2", r"MR\) is such a method", None, None, ""),
    ("10.2-5", "Intro p2", r"cannot easily manufacture", r"less susceptible", None, ""),
    ("10.2-6", "Intro p2", r"edge atlas", r"catalogue|catalog", "DECLINED", "R2 #24 - used 3x, accurate; taste."),
    ("10.2-7", "Intro p3", r"Three questions therefore remain open", r"Three related uncertainties", None, ""),
    ("10.2-8", "Intro p3", r"has adjudicated it", r"unified CKM MR framework", None, ""),
    ("10.2-9", "Intro p4", r"Here we address these questions", r"We therefore constructed", None, ""),
    # ---- 10.3 Results
    ("10.3-1", "Results opening", r"earned that reading", r"benchmark", None, ""),
    ("10.3-2", "Results opening", r"before any staged structure was developed", None, None, ""),
    ("10.3-3", "Results opening (reorder)", None, r"We estimated all 132 non-self directed", "DECLINED",
     "Coverage-first reorder declined (R2 #18 / R3 6.2)."),
    ("10.3-4", "Results heading", r"constrained nulls set the claim", r"not beyond a degree-preserving null", None, ""),
    ("10.3-5", "Results stage 3", r"inserts cleanly", r"no backward CAC edge", None, ""),
    # 10.3-6: bare "as feedback" matched sentences that ALREADY offer index-event as an alternative.
    # The reviewable residue is the DEFAULT label with no alternative in the same clause.
    ("10.3-6", "Results ledger", r"as feedback or an index-event effect", r"ascertainment|liability-scale", "PARTIAL",
     "Prose offers feedback OR index-event, but not pleiotropy/ascertainment as the review asks."),
    ("10.3-7", "Results heading", r"diabetes reaches it predominantly through CAD", r"CAD-independent association with", None, ""),
    ("10.3-8", "Results BMI", r"roughly 80% is direct", r"not mediated through", None, ""),
    ("10.3-9", "Results T2D", r"predominantly, not exclusively, through", r"no detectable CAD-independent|compatible with", None, ""),
    ("10.3-10", "Results subtypes", r"pan-subtype|reaches both subtypes comparably", r"equivalence .*not established|absence of evidence", None, ""),
    ("10.3-11", "Results subtypes", r"HFrEF-restricted", r"P = 2\.2|between-subtype|heterogeneity", None, ""),
    ("10.3-12", "Results heading", r"vicious", r"favours CAD|CAD.{0,3}HF", None, "bare v-cycles regex missed the phrase vicious CAD-HF cycle"),
    ("10.3-13", "Results HF-CAD", r"too thin", r"too imprecise", None, ""),
    ("10.3-14", "Results CAUSE", r"CAUSE backed", r"favoured the causal model|causal model over the sharing", None, ""),
    ("10.3-15", "Results cross-ancestry", r"track the identity line", r"directional(ly)? concordan", None, ""),
    ("10.3-16", "Results TPMI", r"BBJ null overturned|overturned by TPMI", r"aligns with the EUR", None, "bare 'overturn' matched the falsifiability sentence"),
    ("10.3-17", "Results lean diabetes", r"not a selection artifact", r"ascertainment", None, ""),
    ("10.3-18", "Results WHRadjBMI", r"genuine null", r"no detectable WHRadjBMI|no WHRadjBMI.{0,3}HF association", None, ""),
    ("10.3-19", "Results lipids", r"The last ambiguity lay in", None, None, ""),
    # ---- 10.4 Discussion / Conclusions
    ("10.4-1", "Discussion opening", None, r"CAD-independent", "DECLINED",
     "Wholesale opening rewrite; substance already carried."),
    ("10.4-2", "Discussion BMI", r"direct driver", r"beyond CAD|CAD-independent", None, ""),
    ("10.4-3", "Discussion mediation", r"about 80% direct", r"not mediated through", None, ""),
    ("10.4-4", "Discussion subtypes", r"pan-subtype", None, None, ""),
    ("10.4-5", "Discussion T2D", r"coronary-mediated limb", r"without proving|not exclusive|imprecis", None, ""),
    ("10.4-6", "Discussion WHR", r"mass-first story", r"conditioned phenotype|heterogeneity", None, ""),
    ("10.4-7", "Discussion staging", r"internally coherent", r"descriptively concordant", None, ""),
    ("10.4-8", "Discussion portability", r"topology is portable|topology travels", r"suggestive|includes chance", None, ""),
    ("10.4-9", "Discussion overlap", r"directions are not\b", r"magnitude", None, ""),
    ("10.4-10", "Conclusion", r"direct driver of heart failure", r"CAD-independent", None, ""),
    ("10.4-11", "Conclusion staging", r"corroborates the AHA staging", r"no more concordant than degree-matched|degree-preserving", None, ""),
    # Verified by reading CONCLUSIONS.md 2026-07-26: the lean-diabetes paragraph is still there.
    # Deliberately out of scope for v4 (author decision), not an oversight.
    ("10.4-12", "Conclusion lean diabetes", r"lean-diabetes signature", None, None,
     "still in CONCLUSIONS; review asked to move it to the Discussion - open author decision"),
    # ---- 10.5 Figure labels
    ("10.5-1", "Fig 1a", r"resolved in Fig", r"evaluated in Figure 3", None, ""),
    ("10.5-2", "Fig 1c", None, None, "DECLINED", "Facet split of 1c by outcome scale."),
    ("10.5-3", "Fig 1e", r"per-SD / log-OR", None, None, ""),
    ("10.5-4", "Fig 1g", r"inserts cleanly", r"no backward CAC edge", None, ""),
    ("10.5-5", "Fig 2a", r"~80% DIRECT|80% DIRECT", r"not mediated", None, ""),
    ("10.5-6", "Fig 2a", r"predominantly via CAD", r"predominantly via CAD", None, ""),
    ("10.5-7", "Fig 2c", r"fully mediated", r"unstable", None, ""),
    # Fixed 2026-07-26: panel asserted a subtype difference the formal test supports only nominally.
    ("10.5-8", "Fig 2f", r"T2D reaches HFrEF,\s*not HFpEF", r"HFpEF not", None,
     "panel now states the nominal-only caveat; banned string is asserted by gate_package Gate F4"),
    ("10.5-9", "Fig 3b", r"sharing or pleiotropy", r"causal model favoured", None, ""),
    ("10.5-10", "Fig 3e", r"false positive|true positive", r"positive.control|correlated.marker", None, ""),
    ("10.5-11", "Fig 4a", r"not sig \(power\)|power artifact", r"not significant", None, ""),
    ("10.5-12", "Fig 4c", None, r"0\.436|global offset|EAS = 0", "DONE-DIFF",
     "Review wanted sign/z concordance primary. Done differently: fig4.R draws the FITTED offset "
     "(GLOBAL_SLOPE, read from interaction_global_scale.csv) with the identity line de-emphasised."),
    ("10.5-13", "Fig 4e", r"BBJ null overturned|overturned by TPMI", r"aligns with the EUR", None, "same tightening as 10.3-16"),
    ("10.5-14", "Fig 4f", r"[Nn]ode-level bootstrap", r"node-block|resampling", None, ""),
    ("10.5-15", "Graphical abstract", r"DIRECT \+0\.41", r"CAD-independent", None, ""),
    ("10.5-16", "SupplFig7", r"68.{0,6}79\s*%", r"imprecis|unstable", None, "PDF en-dash is \226; old regex false-negatived a REAL stale number"),
    # Fixed 2026-07-26 (fig6.R tier vector). Now asserted continuously by gate_package Gate F3,
    # which requires EXACTLY ONE "HIGHER CONF." row in SupplFig7.pdf - the BMI row.
    ("10.5-17", "SupplFig7", None, r"HIGHER CONF", "DONE",
     "staging demoted to SUPPORTIVE; Gate F3 asserts exactly 1 HIGHER CONF. row in SupplFig7.pdf"),
    ("10.5-18", "SupplFig7", r"not a selection artifact", None, None, ""),
]

# 6.3 diction list - the 17 phrases
DICTION = ["direct driver", "predominantly via", r"\bresolved\b", r"\bcleanly\b", "genuine null",
           "true causal", "false positive", "not a selection artifact", "pan-subtype",
           "track the identity line", "vicious", "load-bearing", "earned that reading",
           "hardened the hedge", "knife-edge"]


def hits(rx):
    return sorted({k for k, t in SURF.items() if re.search(rx, t, re.I)})


print(f"Surfaces loaded: {len(SURF)}\n")
print("=" * 100)
print("SECTION 10 - PROSE REVISION ROWS")
print("=" * 100)
tally = {}
for rid, sec, orig, adopt, override, note in ROWS:
    oh = hits(orig) if orig else []
    ah = hits(adopt) if adopt else []
    if override in ("DECLINED", "CHECK", "PARTIAL", "DONE", "DONE-DIFF"):
        v = override
    elif oh:
        v = "NOT DONE"
    elif adopt and ah:
        v = "DONE"
    elif adopt:
        v = "DONE-DIFF"
    else:
        v = "DONE"
    tally[v] = tally.get(v, 0) + 1
    print(f"[{v:9s}] {rid:8s} {sec}")
    if orig:
        print(f"            original {'STILL PRESENT in ' + str(oh[:5]) if oh else 'absent'}")
    if adopt:
        print(f"            adopted  {'present in ' + str(ah[:4]) if ah else 'NOT FOUND'}")
    if note:
        print(f"            note: {note}")

print("\n" + "=" * 100)
print("SECTION 6.3 - DICTION (phrases the review asks to remove)")
print("=" * 100)
for d in DICTION:
    h = hits(d)
    print(f"[{'HIT ' if h else 'clean'}] {d:32s} {h[:5] if h else ''}")

print("\n" + "=" * 100)
print("TALLY: " + " | ".join(f"{k}={v}" for k, v in sorted(tally.items())))
print("=" * 100)
