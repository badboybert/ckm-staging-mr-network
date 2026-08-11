# -*- coding: utf-8 -*-
"""SupplFig6 data prep: capture the label-permutation NULL distribution of cross-stage
concordance for the EUR / BBJ / TPMI causal networks.

This is a light, non-destructive COPY of the staging logic in scripts/07_staging.py
(same STAGE dict, same Bonferroni-sig + Steiger-correct DAG, same seed 12345, same NP=10000),
modified ONLY to RETURN the null concordance array so it can be histogrammed.
Original 07_staging.py is not modified.

Outputs (figures/suppfig_data/):
  suppfig6_null_<cohort>.csv   -- one column 'concordance' = the 10,000 permuted concordances
  suppfig6_observed.csv        -- cohort, concordance, perm_p, n_cross, n_fwd, n_back, n_tested
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

import csv, io, os, random

BASE = P4_BASE
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(BASE, "figures", "suppfig_data")
os.makedirs(OUT, exist_ok=True)

# AHA 2023 CKM stage per node (identical to 07_staging.py)
STAGE = {"BMI":1, "SBP":2,"TG":2,"HDL":2,"TC":2,"LDL":2,"HbA1c":2,"FI":2,"T2D":2,"eGFR":2,"CKD":2,
         "NAFLD":2, "CAD":4,"HF":4,"Stroke":4,"AF":4}

def load_all(path):
    rows=[]
    with io.open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                p=float(r["ivw_p"]); b=float(r["ivw_b"])
            except Exception:
                continue
            steiger = str(r.get("steiger_correct", r.get("steiger",""))).upper()=="TRUE"
            rows.append((r["exposure"], r["outcome"], b, p, steiger))
    return rows

def analyse(path):
    allrows = load_all(path)
    n_tested = len(allrows)
    bonf = 0.05/n_tested
    # DAG = Bonferroni-significant, Steiger-correct edges among staged nodes
    edges = [(a,b) for (a,b,beta,p,st) in allrows if p<bonf and st and a in STAGE and b in STAGE]
    # de-dup edge set (mirror nx.DiGraph: one edge per ordered pair)
    edge_set = list(dict.fromkeys(edges))
    nodes = list(dict.fromkeys([n for e in edge_set for n in e]))

    cross    = [(a,b) for a,b in edge_set if STAGE[a]!=STAGE[b]]
    forward  = [(a,b) for a,b in cross if STAGE[a]<STAGE[b]]
    backward = [(a,b) for a,b in cross if STAGE[a]>STAGE[b]]
    obs = len(forward)/len(cross) if cross else float('nan')

    # label-permutation null, EXACT: enumerate the distinct stage-label assignments (small support:
    # 1,320 for BBJ's 11 nodes, 840 for TPMI's 9) and take the exact tail, matching 07_staging.py.
    # The plotted null is the exact distribution, not a Monte-Carlo sample of it.
    from collections import Counter as _C
    def _msp(c):
        n=sum(c.values()); cur=[None]*n
        def rec(i):
            if i==n:
                yield tuple(cur); return
            for v in sorted(c):
                if c[v]:
                    c[v]-=1; cur[i]=v
                    yield from rec(i+1)
                    c[v]+=1
        return rec(0)
    null=[]
    for pm in _msp(dict(_C(STAGE[n] for n in nodes))):
        smap=dict(zip(nodes,pm))
        cc=[(a,b) for a,b in edge_set if smap[a]!=smap[b]]
        ff=sum(1 for a,b in cc if smap[a]<smap[b])
        null.append(ff/len(cc) if cc else 0)
    ge=sum(1 for c in null if c>=obs-1e-12)
    perm_p=ge/len(null)                              # EXACT tail (no add-one)
    return dict(concordance=obs, perm_p=perm_p, null=null, n_cross=len(cross),
                n_fwd=len(forward), n_back=len(backward), n_tested=n_tested,
                backward=backward)

COHORTS = [
    ("EUR",  os.path.join(RES, "forward_local_edges.csv")),
    ("BBJ",  os.path.join(RES, "network_eas_edges.csv")),
    ("TPMI", os.path.join(RES, "network_tpmi_full.csv")),
]

summary=[]
for name, path in COHORTS:
    r = analyse(path)
    # write null array
    with io.open(os.path.join(OUT, f"suppfig6_null_{name}.csv"), "w", newline="", encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["concordance"])
        for c in r["null"]: w.writerow([c])
    summary.append((name, r["concordance"], r["perm_p"], r["n_cross"], r["n_fwd"], r["n_back"], r["n_tested"]))
    print(f"{name:5s}  concordance={r['concordance']:.3f}  perm_p={r['perm_p']:.4f}  "
          f"cross={r['n_cross']} (fwd {r['n_fwd']}/back {r['n_back']})  n_tested={r['n_tested']}  "
          f"backward={r['backward']}")

with io.open(os.path.join(OUT, "suppfig6_observed.csv"), "w", newline="", encoding="utf-8") as f:
    w=csv.writer(f); w.writerow(["cohort","concordance","perm_p","n_cross","n_fwd","n_back","n_tested"])
    for row in summary: w.writerow(row)
print("\nWrote suppfig6_observed.csv + 3 null arrays to", OUT)
