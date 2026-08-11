# -*- coding: utf-8 -*-
"""Deterministic guard: confirm the tightened section preserved EVERY numeric token, figure callout,
and citation vs the original (multiset equality). Independent of the verifier agents (verify-the-verifier)."""
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

import re, io, sys
from collections import Counter
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
FIG = os.path.join(P4_BASE, "figures")

SUP = "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"
NUM = re.compile(r"[+\-−]?\d[\d.,]*(?:\s*[×x]\s*10[" + SUP + r"]+)?|10[" + SUP + r"]+|[" + SUP + r"]{2,}")
CALL = re.compile(r"(?:Figure|Fig\.|Supplementary Figure|Supplementary Table)\s*S?\d+[a-f]?(?:[,–\-][a-f])*")
CITE = re.compile(r"[A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+ et al\.? \(?\d{4}\)?|[A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+ and [A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+ \d{4}")

def toks(path, rx):
    t = open(path, encoding="utf-8").read()
    t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
    return Counter(m.group(0).strip() for m in rx.finditer(t))

def norm_num(c):
    # normalise whitespace inside "× 10" (so "1 × 10⁻³" == "1 ×10⁻³") and strip glued edge punctuation
    # (a mid-sentence "3," becoming an end-of-sentence "3." after a split is the SAME value 3).
    return Counter(re.sub(r"\s+", "", k).strip(".,") for k in c.elements())

import os
ok = True
# All prose surfaces, not three of them. A guard whose scope lags the manuscript proves nothing
# about the sections it never opened -- and the 2026-07-25 fluency pass touched CONCLUSIONS and
# FIGURE_LEGENDS, neither of which was in the original list.
for sec in ["INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION", "CONCLUSIONS", "FIGURE_LEGENDS"]:
    orig, tight = f"{FIG}/{sec}.md", f"{FIG}/{sec}.tightened.md"
    if not os.path.exists(tight):
        continue   # already applied + .tightened removed; nothing to compare
    print(f"\n===== {sec} =====")
    for label, rx, nrm in [("NUMBERS", NUM, True), ("CALLOUTS", CALL, False), ("CITATIONS", CITE, False)]:
        a = toks(orig, rx); b = toks(tight, rx)
        if nrm: a, b = norm_num(a), norm_num(b)
        dropped = a - b        # in original, missing/short in tightened
        added   = b - a        # new in tightened
        status = "OK" if not dropped and not added else "MISMATCH"
        if dropped or added: ok = False
        print(f"  {label:9s}: {sum(a.values())} orig / {sum(b.values())} tightened -> {status}")
        if dropped: print(f"      DROPPED (in orig, not in tightened): {dict(dropped)}")
        if added:   print(f"      ADDED   (in tightened, not in orig): {dict(added)}")

print("\n" + ("ALL CLEAN — safe to apply." if ok else "*** MISMATCH — DO NOT APPLY until resolved. ***"))
