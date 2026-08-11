# -*- coding: utf-8 -*-
"""Objective prose scan: sentence lengths, semicolon density (main text) + legend-title lengths."""
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
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
FIG = os.path.join(P4_BASE, "figures")
ABBR = ["et al.", "e.g.", "i.e.", "vs.", "Fig.", "cf.", "No.", "approx.", "ref.", "Dr."]

def sentences(txt):
    t = re.sub(r"<!--.*?-->", "", txt, flags=re.S)
    t = re.sub(r"^#.*$", "", t, flags=re.M)
    t = " ".join(t.split())
    for a in ABBR:
        t = t.replace(a, a.replace(".", ""))
    parts = re.split(r"(?<=[.;])\s+(?=[A-Z(–])", t)
    return [p.replace("", ".").strip() for p in parts if p.strip()]

print("========== MAIN-TEXT SENTENCE LENGTHS ==========")
for name in ["INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION"]:
    txt = open(f"{FIG}/{name}.md", encoding="utf-8").read()
    ss = sentences(txt)
    wl = [(len(s.split()), s) for s in ss]
    longs = sorted([x for x in wl if x[0] >= 38], reverse=True)
    semis = sum(s.count(";") for s in ss)
    mean = sum(w for w, _ in wl) // max(1, len(wl))
    print(f"\n=== {name}: {len(ss)} sentences | mean {mean}w | {len(longs)} long (>=38w) | {semis} semicolons ===")
    for w, s in longs[:7]:
        print(f"  [{w}w] {s[:175]}")

print("\n\n========== FIGURE-LEGEND TITLES ==========")
for fn, tag in [("FIGURE_LEGENDS", "Figure"), ("SUPPL_FIGURES", "Supplementary Figure")]:
    txt = open(f"{FIG}/{fn}.md", encoding="utf-8").read()
    for m in re.finditer(r"\*\*" + tag + r" S?\d[^|]*\|\s*(.+?)\.?\*\*", txt):
        title = m.group(1).strip().rstrip(".")
        w = len(title.split())
        flags = []
        if ";" in title: flags.append("SEMICOLON")
        if ":" in title: flags.append("COLON")
        if w > 12: flags.append(f"{w}w")
        mark = ("  <-- " + ", ".join(flags)) if flags else ""
        head = re.search(tag + r" (S?\d)", m.group(0)).group(1)
        print(f"  {tag[:9]} {head}: \"{title}\"{mark}")
