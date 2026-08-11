# -*- coding: utf-8 -*-
"""Deterministic prose-quality + AI-tell scan over the shipped manuscript sections.

Two families of signal, both computed, neither a substitute for a human/agent read:
  PROSE   sentence-length distribution, long-sentence offenders, passive-voice rate, hedge density.
  AI-TELL formulaic AI vocabulary, "not only ... but also" scaffolding, em-dash density, uniform
          sentence length (low burstiness), and the tricolon "X, Y, and Z" cadence LLMs overuse.

Reports per section; flags offenders. Run: python scripts/prose_ai_scan.py
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

import io, os, re, sys, statistics
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# The package version is DERIVED from build_submission_v1.py, never typed here. This file was
# pinned to "submission package v2" and stayed pinned through two later cuts, so running it would
# have silently scanned a two-versions-stale artifact and reported it as the manuscript. A prose
# gate that reads the wrong package is worse than no prose gate: it manufactures confidence.
_ROOT = P4_ROOT
_bsv = io.open(os.path.join(_ROOT, "independent_build", "scripts", "build_submission_v1.py"),
               encoding="utf-8").read()
_m = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _m, "could not read PKG_VERSION from build_submission_v1.py"
PKG = os.path.join(_ROOT, f"submission package {_m.group(1)}")
assert os.path.isdir(PKG), f"package not found: {PKG}"

# Every prose surface an editor or reviewer actually reads. CONCLUSIONS and the two CVD-required
# additions were absent from this list because they did not exist when it was written; a scan whose
# scope lags the manuscript reports a clean bill for text it never opened.
SECTIONS = ["ABSTRACT", "INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION", "CONCLUSIONS",
            "FIGURE_LEGENDS"]

# AI-tell vocabulary: words/phrases that appear at anomalous rates in LLM prose.
AI_WORDS = [r"\bdelve\b", r"\bdelved\b", r"\btapestry\b", r"\bunderscore[sd]?\b", r"\bunderscoring\b",
            r"\bshowcas(e|es|ed|ing)\b", r"\bleverag(e|es|ed|ing)\b", r"\bnavigat(e|es|ed|ing)\b",
            r"\brealm\b", r"\bintricate\b", r"\bmyriad\b", r"\bpivotal\b", r"\bcrucial\b",
            r"\bnuanced\b", r"\bholistic\b", r"\bmultifaceted\b", r"\bparadigm\b",
            r"\btestament\b", r"\bin the realm of\b", r"\bit is (important|worth noting) to\b",
            r"\bplays? a (crucial|vital|pivotal|key|significant) role\b", r"\brich (tapestry|landscape)\b",
            r"\bever-(evolving|changing|growing)\b", r"\bcornerstone\b", r"\bharness(e|es|ed|ing)?\b",
            r"\bunveil(s|ed|ing)?\b", r"\bshed(s|ding)? light\b", r"\bat the intersection of\b",
            r"\bnavigat", r"\bfoster(s|ed|ing)?\b", r"\bbolster(s|ed|ing)?\b", r"\brobustly\b"]
SCAFFOLD = [r"\bnot only\b[^.]{0,80}\bbut (also|rather)\b", r"\bit is worth noting that\b",
            r"\bit is important to (note|recognize|emphasize)\b", r"\bin conclusion\b",
            r"\bfurthermore\b", r"\bmoreover\b", r"\badditionally\b", r"\bin summary\b",
            r"\bfirstly\b", r"\bsecondly\b", r"\bthat being said\b"]


def sentences(text):
    text = re.sub(r"\s+", " ", text)
    # protect common abbreviations and "P = ..." from being split
    text = text.replace("et al.", "et al<DOT>").replace("e.g.", "e<DOT>g<DOT>").replace("i.e.", "i<DOT>e<DOT>")
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(])", text)
    return [p.replace("<DOT>", ".").strip() for p in parts if len(p.strip()) > 3]


def wordcount(s):
    return len(re.findall(r"[A-Za-z][A-Za-z'-]*", s))


def read_section(name):
    t = io.open(os.path.join(PKG, "01_manuscript", name + ".md"), encoding="utf-8").read()
    t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
    t = "\n".join(ln for ln in t.split("\n")
                  if not ln.startswith("#") and not ln.strip().startswith("**Title")
                  and not ln.strip() in ("---", "") and not ln.strip().startswith("## "))
    # strip markdown emphasis and inline math markers
    t = re.sub(r"[*_`]", "", t)
    return t


print("=" * 96)
print(f"PROSE + AI-TELL SCAN — {os.path.basename(PKG)}, {len(SECTIONS)} shipped sections")
print("=" * 96)
all_lens = []
for name in SECTIONS:
    t = read_section(name)
    ss = sentences(t)
    lens = [wordcount(s) for s in ss]
    if not lens:
        continue
    all_lens += lens
    mean = statistics.mean(lens)
    sd = statistics.pstdev(lens)
    long40 = [s for s in ss if wordcount(s) > 40]
    passive = sum(1 for s in ss if re.search(r"\b(was|were|is|are|been|be|being)\s+\w+(ed|en)\b", s))
    ai_hits = [(re.search(p, t, re.I).group(0), len(re.findall(p, t, re.I)))
               for p in AI_WORDS if re.search(p, t, re.I)]
    scaf = [(re.search(p, t, re.I).group(0)[:40], len(re.findall(p, t, re.I)))
            for p in SCAFFOLD if re.search(p, t, re.I)]
    emdash = t.count("—")
    tricolon = len(re.findall(r"\w+, \w[\w -]*, and \w", t))
    print("\n### %-13s  %d sentences | mean %.1f w | SD %.1f (burstiness) | %d words"
          % (name, len(ss), mean, sd, sum(lens)))
    print("    long (>40w): %d   passive-ish: %d/%d (%.0f%%)   em-dashes: %d   tricolons: %d"
          % (len(long40), passive, len(ss), 100 * passive / len(ss), emdash, tricolon))
    if ai_hits:
        print("    AI-tell words: " + ", ".join("%s(%d)" % (w, n) for w, n in ai_hits))
    else:
        print("    AI-tell words: none")
    if scaf:
        print("    scaffolding:   " + ", ".join("%s(%d)" % (w, n) for w, n in scaf))
    for s in sorted(long40, key=wordcount, reverse=True)[:3]:
        print("      LONG(%dw): %s..." % (wordcount(s), s[:150]))

print("\n" + "=" * 96)
gm = statistics.mean(all_lens); gs = statistics.pstdev(all_lens)
print("WHOLE MANUSCRIPT: %d sentences | mean %.1f w | SD %.1f | CV %.2f"
      % (len(all_lens), gm, gs, gs / gm))
print("Reference bands: human scientific prose mean ~18-25 w, SD ~10-14 (CV ~0.5-0.65).")
print("Low burstiness (CV < 0.4) is the single strongest statistical AI-tell; high CV is human-like.")
