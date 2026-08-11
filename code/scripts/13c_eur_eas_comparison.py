# -*- coding: utf-8 -*-
"""Step 13c: regenerate results/eur_vs_eas_comparison.csv, the EUR-vs-EAS edge comparison table.

This file is read by 14_standardise.py, figures/fig4.R, figures/prep_fig5.py and the supplementary
workbook, and it ships as source data — but no script in the tree wrote it, so its provenance could
not be checked. This reconstructs it from the two committed network tables and asserts that the
result matches the committed file, so the artifact stops being hand-held.

Columns are unchanged: edge, EUR_b, EUR_p, EAS_b, EAS_p, dir_concordant, both_sig.
Betas and P-values are written at the same 3-significant-figure precision as the committed file;
downstream consumers that need full precision take it from the network tables directly.
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

import csv, io, os, sys

BASE = P4_BASE
RES = os.path.join(BASE, "results")
OUT = os.path.join(RES, "eur_vs_eas_comparison.csv")
EUR_BONF_N = None          # derived below (reported, not used as the both_sig rule)


def load(path):
    d = {}
    with io.open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                d[(r["exposure"], r["outcome"])] = (float(r["ivw_b"]), float(r["ivw_p"]))
            except (KeyError, ValueError):
                continue
    return d


eur = load(os.path.join(RES, "forward_local_edges.csv"))
eas = load(os.path.join(RES, "network_eas_edges.csv"))
EUR_BONF_N = len(eur)
BONF = 0.05 / EUR_BONF_N

shared = [k for k in eur if k in eas]
rows = []
for a, b in sorted(shared):
    eb, ep = eur[(a, b)]
    ab, ap = eas[(a, b)]
    rows.append(dict(
        edge="%s->%s" % (a, b),
        EUR_b="%.3f" % eb, EUR_p="%.1e" % ep,
        EAS_b="%.3f" % ab, EAS_p="%.1e" % ap,
        dir_concordant=str((eb > 0) == (ab > 0)),
        # both_sig marks edges nominally significant in BOTH ancestries (P < 0.05 each side); it is a
        # display flag for Figure 4c, not the Bonferroni-significant subset used for portability.
        both_sig=str(ep < 0.05 and ap < 0.05)))

# --- reproduce-the-committed-file check ---------------------------------------------------------
if os.path.exists(OUT):
    old = {r["edge"]: r for r in csv.DictReader(io.open(OUT, encoding="utf-8"))}
    new = {r["edge"]: r for r in rows}
    only_old = sorted(set(old) - set(new))
    only_new = sorted(set(new) - set(old))
    diff = [(e, old[e], new[e]) for e in sorted(set(old) & set(new))
            if any(abs(float(old[e][c]) - float(new[e][c])) > 5e-3 for c in ("EUR_b", "EAS_b"))
            or old[e]["dir_concordant"] != new[e]["dir_concordant"]
            or old[e]["both_sig"] != new[e]["both_sig"]]
    print("committed rows %d | regenerated rows %d" % (len(old), len(new)))
    print("edges only in the committed file: %s" % (only_old or "none"))
    print("edges only in the regenerated file: %s" % (only_new or "none"))
    print("rows whose beta/verdict differ: %d" % len(diff))
    for e, o, n in diff[:10]:
        print("   %-14s committed %s/%s/%s  regenerated %s/%s/%s"
              % (e, o["EUR_b"], o["EAS_b"], o["both_sig"], n["EUR_b"], n["EAS_b"], n["both_sig"]))
    if only_old or only_new or diff:
        print("\nREFUSING TO OVERWRITE: the regeneration does not reproduce the committed file.")
        print("Investigate before adopting — pass --force to overwrite deliberately.")
        if "--force" not in sys.argv:
            sys.exit(1)

with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["edge", "EUR_b", "EUR_p", "EAS_b", "EAS_p",
                                       "dir_concordant", "both_sig"])
    w.writeheader()
    for r in rows:
        w.writerow(r)
print("\nwrote %s (%d shared edges; EUR Bonferroni 0.05/%d)" % (OUT, len(rows), EUR_BONF_N))
