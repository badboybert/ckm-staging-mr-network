# -*- coding: utf-8 -*-
"""TRUTH-ORIENTED claim scan over every authored source (not a regression list).

Three families, all of which shipped in package v2 and none of which qa_manuscript.py or
gate_package.py Gate D0 could see:

  R  REVISION HISTORY — this is a first submission. No shipped file may narrate an earlier version,
     an internal review, or a rebuild ("in response to peer review", "we therefore rebuilt", "the
     earlier formulation", reviewer item codes such as C4/E3).
  C  RETRACTED CLAIMS — wording the manuscript's own recalibrations removed, which survived in the
     cover letter, the STROBE checklist and the workbook because no gate read them.
  S  SUPERSEDED NUMBERS — values from the 111-edge network or the retired ledger.

Run:  python scripts/scan_claims.py            (scans the authored sources in figures/)
      python scripts/scan_claims.py <dir>      (scans any tree, e.g. a cut package)
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

import io, os, re, sys, glob

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = P4_BASE

PATTERNS = [
    # ---- R: revision history / internal-review provenance -------------------------------------
    ("R", "in response to peer review", r"in response to (peer )?review"),
    ("R", "narrates a rebuild", r"we (therefore )?rebuilt|was rebuilt|has been rebuilt"),
    ("R", "narrates an earlier version", r"[Tt]he earlier formulation|the previous (version|draft|formulation)|in the earlier (draft|version)"),
    ("R", "supersession language", r"\bsupersed(e|es|ed|ing)\b"),
    ("R", "reviewer item code", r"\(reviewer [CEMH]\d"),
    # Added 2026-07-25. The internal review is INTERNAL: no shipped file may name it, cite it or
    # attribute a change to it. This escaped every earlier rule because the phrase was "reviewer
    # round-2 showed", inside a WORKBOOK caption -- a surface no prose gate had ever read. The
    # ban is on the actor, not on a particular sentence, so paraphrases cannot slip through.
    ("R", "names the internal review", r"\breviewers?\b(?!\s|$)|\breviewer (round|comment|item|report)|round[- ]?[12] (review|response)|rebuttal|response to (the )?review"),
    ("R", "revision round", r"\b(round[- ]1|first round|this revision|during revision)\b"),
    # ---- C: retracted claims -------------------------------------------------------------------
    ("C", "reproduces the staging order", r"reproduce[sd]? the AHA staging order|reproduces the order"),
    ("C", "topology ports/replicates", r"topology (replicate|port)|topology (ports|replicates)|topology carr(y|ies|ied) to"),
    ("C", "replicated ancestry divergence", r"replicated .{0,30}ancestry divergence|replicated systolic|divergence we replicate"),
    ("C", "does not port", r"does not port"),
    ("C", "genuine feedback", r"genuine feedback|genuine, not artifactual|genuine reciprocal"),
    ("C", "licenses prevention", r"\blicens(e|es|ed)\b"),
    ("C", "only via CAD", r"only via CAD|reaches heart failure only|only through coronary"),
    ("C", "fully mediated", r"fully mediated"),
    ("C", "pre-registered", r"[Pp]re-?registered"),
    ("C", "AHA-2023 as a live construct", r"AHA-2023"),
    ("C", "full battery on every edge", r"[Ee]very (headline )?edge (in the .{0,20}network )?(was|carries|carried) .{0,40}(full|fixed) .{0,20}battery"),
    ("C", "retired title", r"Adiposity is a common driver of the cardiovascular"),
    # Retired 2026-07-25 (v3.2) when the title moved to the adiposity-centred architecture wording.
    # Scoped to the distinctive "map of selected" clause: the phrase "adiposity as a direct driver of
    # heart failure" is still LEGITIMATE prose in the abstract's Conclusions and the Results synthesis,
    # so banning the whole old title string would be wrong.
    ("C", "retired title (v3 map-of-selected form)", r"Mendelian-randomization map of selected"),
    # Added 2026-07-25 after checking the round-2 review's wording demands against the SHIPPED v3.2
    # package rather than against the fix list. Three classes had survived the earlier sweeps:
    #  - "negative control" without the correlated-marker qualifier. HDL->CAD is not a clean null;
    #    the qualifier is what carries that, so a bare mention undoes the concession.
    ("C", "unqualified negative control", r"(?<!correlated-marker )negative control"),
    #  - the "topology travels" slogan, which survived in a Results section heading.
    ("C", "topology-travels slogan", r"topology travels"),
    #  - the ASSERT/DEVELOP/BOUND tier vocabulary in PROSE. The panel labels were changed to
    #    HIGHER CONF./SUPPORTIVE/EXPLORATORY but the text kept the verbs, so figure and text
    #    disagreed -- the exact cross-file class the review's consistency audit is about.
    ("C", "retired evidence-tier vocabulary", r"evidence asserts|develops as hypothesis|asserts, develops"),
    # ---- S: superseded numbers ------------------------------------------------------------------
    ("S", "0.05/111", r"0\.05/111|4\.5 ?× ?10⁻⁴|4\.5e-0?4"),
    ("S", "111 as the current network", r"111[- ]edge|111 directed|\bthe 111\b"),
    # Round-7: scoped to a STAGING context. The bare literal also matches a legitimate edge-level
    # P — the S5 legend reports HF→T2D at P = 1.8 × 10⁻³ against the Mahajan-noUKB outcome — and a
    # rule that fires on any occurrence of four digits stops being about staging at all.
    ("S", "retired staging P (1.8e-3)",
     r"(?:permutation|concordance|staging|0\.926)[^.]{0,80}(?:1\.8 ?× ?10⁻³|(?<![\d.])0\.0018(?![\d]))"
     r"|(?:1\.8 ?× ?10⁻³|(?<![\d.])0\.0018(?![\d]))[^.]{0,60}(?:permutation|staging|concordance)"),
    ("S", "superseded MC staging P (1.2e-3, now exact 1.5e-3)", r"(?:permutation|concordance|0\.926|staging)[^.]{0,60}1\.2 ?× ?10⁻³|1\.2 ?× ?10⁻³[^.]{0,40}(?:permutation|Figure 1b)"),
    ("S", "retired ledger counts", r"11 concordant|eleven concordant"),
    ("S", "Kendall tau portability", r"Kendall|τ ?= ?0\.49|44 ?(of|/) ?65"),
    ("S", "noUKB staging 1.000/1.1e-3", r"1\.000, ?\*?P\*? ?= ?1\.1 ?× ?10⁻³|0\.957, ?\*?P\*? ?= ?3 ?× ?10⁻⁴"),
    # Added 2026-07-24, each for a defect that shipped past every existing rule:
    ("S", "noUKB staging Monte-Carlo P (3e-4; not attainable on the exact 1,980-labeling support)",
     r"0\.960[^.]{0,40}3 ?× ?10⁻⁴|0\.913[^.]{0,40}3 ?× ?10⁻⁴|0\.9(?:60|13), ?P ?= ?0\.0003"),
    ("S", "pre-enumeration CAC staging P (1e-4, now exact 1.2e-4)",
     r"0\.941[^.]{0,40}(?<!\.)1 ?× ?10⁻⁴"),
    # A retired METHOD description, not a number. The staging null is now an exhaustive enumeration;
    # the 10,000-replicate bootstrap and sign-flip null are separate and legitimate, so this rule is
    # scoped to the label-permutation/staging test only.
    ("C", "staging null described as Monte-Carlo",
     r"(?:label-permutation|staging)[^.]{0,60}10,?000|10,?000[- ](?:shuffle|iteration)[^.]{0,40}(?:label-permutation|staging)"),
]

# Legitimate historical statements, matched verbatim and exempted with a reason.
ALLOW = [
    ("survived completing the network from 111 to 132 edges",
     "true historical statement about network completion, not a live claim"),
    ("from 111 to 132",
     "true historical statement about network completion"),
    ("*P* = 4.7 × 10⁻⁴ and 1.8 × 10⁻³",
     "these are the CAD->T2D and HF->T2D noUKB edge P-values, not the staging permutation P"),
    ("is scanned for superseded numbers, retracted wording and revision-history language",
     "the README describing what this very gate checks — a statement ABOUT the scan, not a claim in the package"),
]


def read_any(p):
    """Text surface of a shipped file. .docx and .xlsx are NOT optional: the cover letter, the STROBE
    checklist and the supplementary workbook are where retracted claims survived a .md-only scan."""
    e = os.path.splitext(p)[1].lower()
    if e in (".md", ".txt", ".csv", ".json"):
        try:
            return io.open(p, encoding="utf-8", errors="replace").read()
        except (IsADirectoryError, OSError):
            return None
    if e == ".docx":
        from docx import Document
        d = Document(p)
        parts = [q.text for q in d.paragraphs]
        for tb in d.tables:
            for row in tb.rows:
                parts += [c.text for c in row.cells]
        return chr(10).join(parts)
    if e == ".xlsx":
        import openpyxl
        wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
        parts = list(wb.sheetnames)
        for ws in wb.worksheets:
            for row in ws.iter_rows(values_only=True):
                parts += [v for v in row if isinstance(v, str)]
        wb.close()
        return chr(10).join(parts)
    return None                                   # png/pdf — handled by the gate's figure scan


def scan(paths):
    fails = 0
    for p in sorted(paths):
        t = read_any(p)
        if t is None:
            continue
        t = re.sub(r"<!--.*?-->", "", t, flags=re.S)          # build comments never ship
        for fam, label, pat in PATTERNS:
            for m in re.finditer(pat, t):
                frag = " ".join(t[max(0, m.start() - 70):m.end() + 70].split())
                # Compare against the allowlist with markdown emphasis stripped: the same sentence
                # appears as "*P* = ..." in the .md and "P = ..." in the rendered .docx, and an
                # allowlist that only matched one form let the .docx fire a known false positive.
                bare = frag.replace("*", "").replace("_", "").replace("`", "")
                if any(a in frag or a.replace("*", "") in bare for a, _ in ALLOW):
                    continue
                fails += 1
                print("  [%s] %-34s %-28s ...%s..." % (fam, label, os.path.basename(p), frag))
    return fails


# The authored sources of the files that actually SHIP. Listed explicitly rather than globbed:
# figures/ also holds superseded drafts and verification notes, and a glob would drown the real
# hits in noise from files no reader ever sees.
SHIPPED_SOURCES = ["ABSTRACT.md", "INTRODUCTION.md", "METHODS.md", "RESULTS.md", "DISCUSSION.md",
                   "FIGURE_LEGENDS.md", "SUPPL_FIGURES.md", "FRONT_MATTER.md", "COVER_LETTER.md",
                   "STROBE_MR_CHECKLIST.md", "REFERENCES_NUMBERED.md"]

if __name__ == "__main__":
    if len(sys.argv) > 1:                       # scan a cut package: every text file in it
        root = sys.argv[1]
        files = []
        for _ext in ("*.md", "*.txt", "*.docx", "*.xlsx"):
            files += glob.glob(os.path.join(root, "**", _ext), recursive=True)
    else:                                       # scan the authored sources of the shipped files
        root = os.path.join(BASE, "figures")
        files = [os.path.join(root, f) for f in SHIPPED_SOURCES]
        missing = [f for f in files if not os.path.exists(f)]
        if missing:
            print("SCAN ABORT: shipped source(s) missing:", missing); sys.exit(2)
    print("scan_claims: %d file(s) under %s" % (len(files), root))
    n = scan(files)
    print("\n===== CLAIM SCAN %s — %d hit(s) =====" % ("CLEAN" if n == 0 else "FAILED", n))
    sys.exit(1 if n else 0)
