# -*- coding: utf-8 -*-
"""Re-derive every package count quoted in the internal response letters from CANONICAL_FACTS.

WHY THIS EXISTS. The letters quote the package's headline counts. Every edit that moves a word
count makes them stale, and chasing the literals by hand has now gone wrong twice in two different
ways: a bulk regex ate the commas separating clauses ("Methods 3,975 Results 4,291"), and targeted
replacements kept missing by one because the literal had already moved since the last run.

So: match the SHAPE of each statement and rewrite its numbers from the manifest. Never match a
specific value -- that is what kept going stale.

TRAP THIS MUST NOT REPEAT: R4 carries "and **54**" as the GRAPH SIZE, not a reference count. Nothing
here touches a bare number; every pattern is anchored on its surrounding words.

Run after derive_facts.py, then rebuild the letter .docx twins.
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

import glob
import io
import json
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = json.load(io.open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))
w = F["word_counts"]
MAIN, REFS, AB = F["main_text_words"], F["references_cited"], F["abstract_full_words"]
BASELINE = 12436                      # the pre-trim main-text length the letters compare against
TRIM = (BASELINE - MAIN) / BASELINE * 100

SUBS = [
    # "main text **11,518** words"  /  "Main text **11,518 words**"
    (re.compile(r"(main text \*\*)[\d,]+(\*\* words)", re.I), lambda m: f"{m.group(1)}{MAIN:,}{m.group(2)}"),
    (re.compile(r"(Main text \*\*)[\d,]+( words\*\*)"), lambda m: f"{m.group(1)}{MAIN:,}{m.group(2)}"),
    # "(Introduction 703, Methods 4,024, Results 4,294, Discussion 2,283, Conclusions 263)"
    (re.compile(r"Introduction [\d,]+, Methods [\d,]+, Results [\d,]+, Discussion [\d,]+, "
                r"Conclusions [\d,]+"),
     lambda m: (f"Introduction {w['INTRODUCTION']:,}, Methods {w['METHODS']:,}, "
                f"Results {w['RESULTS']:,}, Discussion {w['DISCUSSION']:,}, "
                f"Conclusions {w['CONCLUSIONS']:,}")),
    # "12,436 → 11,518 words, −7.0%"
    (re.compile(r"(12,436 → )[\d,]+( words, )−?[\d.]+%"),
     lambda m: f"{m.group(1)}{MAIN:,}{m.group(2)}−{TRIM:.1f}%"),
    # "Abstract **324** words"
    (re.compile(r"(Abstract \*\*)\d{3}(\*\* words)"), lambda m: f"{m.group(1)}{AB}{m.group(2)}"),
    # reference counts, always anchored on the word "references"/"entries ... resolve"
    (re.compile(r"\*\*\d{2,3} references\*\*"), lambda m: f"**{REFS} references**"),
    (re.compile(r"\b\d{2,3}( references in numbered Vancouver order)"),
     lambda m: f"{REFS}{m.group(1)}"),
    (re.compile(r"All \d{2,3}( entries now resolve)"), lambda m: f"All {REFS}{m.group(1)}"),
    (re.compile(r"(references resolved against PubMed/Crossref \| )\d{2,3} / \d{2,3}"),
     lambda m: f"{m.group(1)}{REFS} / {REFS}"),
]


def main():
    print(f"canonical: main {MAIN:,} | abstract {AB} | refs {REFS} | trim −{TRIM:.1f}% | "
          f"sections {dict(w)}")
    total = 0
    for p in sorted(glob.glob(os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R*.md"))):
        t0 = io.open(p, encoding="utf-8").read()
        t, n = t0, 0
        for rx, rep in SUBS:
            t, k = rx.subn(rep, t)
            n += k
        if t != t0:
            io.open(p, "w", encoding="utf-8", newline="\n").write(t)
            print(f"  {os.path.basename(p):38s} {n} rewritten")
            total += n
    print(f"{total} statement(s) re-derived")

    # the graph size must survive untouched -- it is a 54 that is NOT a reference count
    r4 = os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R4.md")
    if os.path.exists(r4):
        assert "and **54**" in io.open(r4, encoding="utf-8").read(), \
            "R4's graph-size 54 was destroyed -- it is not a reference count"
        print("checked: R4's graph-size 54 is intact")

    # and nothing may be left mangled the way the old bulk regex left it
    bad = []
    for p in sorted(glob.glob(os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R*.md"))):
        for i, ln in enumerate(io.open(p, encoding="utf-8").read().split("\n"), 1):
            # NOTE the [\d,]*\d: a greedy [\d,]+ swallows the comma that legitimately separates the
            # clauses, so "Results 4,294, Discussion" matched and the detector cried wolf on
            # correctly-formatted lines. Require the number to END before the separator.
            if re.search(r"Methods [\d,]*\d Results|\*\*\d+\*\* words\*\* words|"
                         r"Results [\d,]*\d Discussion", ln):
                bad.append(f"{os.path.basename(p)}:{i}")
    print(f"mangled lines: {bad or 'none'}")


if __name__ == "__main__":
    main()
