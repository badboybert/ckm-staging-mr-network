# -*- coding: utf-8 -*-
"""Figure 1 data prep: significant EUR causal edges, cross-stage classification,
permutation null for staging concordance, positive controls. -> figures/data/fig1_*.csv"""
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

import csv, io, os, math, random
BASE = P4_BASE
OUT = os.path.join(BASE, "figures/data"); os.makedirs(OUT, exist_ok=True)
STAGE = {"BMI":1,"SBP":2,"TG":2,"HDL":2,"TC":2,"LDL":2,"HbA1c":2,"T2D":2,"eGFR":2,"CKD":2,"CAD":4,"HF":4,"Stroke":4,"AF":4}

rows = []
with io.open(os.path.join(BASE,"results/forward_local_edges.csv"),encoding="utf-8") as f:
    for r in csv.DictReader(f):
        try: b=float(r["ivw_b"]); p=float(r["ivw_p"]); se=float(r["ivw_se"])
        except: continue
        st = str(r.get("steiger_correct","")).upper()=="TRUE"
        ei = r.get("egger_intercept_p",""); qp = r.get("Q_p","")
        rows.append(dict(exp=r["exposure"], out=r["outcome"], b=b, se=se, p=p, st=st,
                         egg_ip=float(ei) if ei not in ("","NA") else float("nan"),
                         qp=float(qp) if qp not in ("","NA") else float("nan"),
                         nsnp=int(float(r["nsnp"]))))
n_tested = len(rows); bonf = 0.05/n_tested

# significant Steiger-correct edges among staged nodes
sig = [r for r in rows if r["p"]<bonf and r["st"] and r["exp"] in STAGE and r["out"] in STAGE]
with io.open(os.path.join(OUT,"fig1_edges.csv"),"w",encoding="utf-8",newline="") as f:
    w=csv.writer(f); w.writerow(["exp","out","b","p","nsnp","stage_from","stage_to","cross","direction","egger_int_p","Q_p"])
    for r in sig:
        sf,st_=STAGE[r["exp"]],STAGE[r["out"]]
        cross = sf!=st_
        direction = "forward" if sf<st_ else ("backward" if sf>st_ else "intra")
        w.writerow([r["exp"],r["out"],f"{r['b']:.4f}",f"{r['p']:.3e}",r["nsnp"],sf,st_,cross,direction,
                    f"{r['egg_ip']:.4f}" if r['egg_ip']==r['egg_ip'] else "NA",
                    f"{r['qp']:.2e}" if r['qp']==r['qp'] else "NA"])
cross=[r for r in sig if STAGE[r["exp"]]!=STAGE[r["out"]]]
fwd=[r for r in cross if STAGE[r["exp"]]<STAGE[r["out"]]]
back=[r for r in cross if STAGE[r["exp"]]>STAGE[r["out"]]]
conc = len(fwd)/len(cross)
print(f"n_tested={n_tested} bonf={bonf:.2e} | sig-staged edges={len(sig)} | cross={len(cross)} fwd={len(fwd)} back={len(back)} concordance={conc:.3f}")
for r in back:
    print(f"  backward {r['exp']}->{r['out']} b={r['b']:.3f} eggIp={r['egg_ip']:.3f} Qp={r['qp']:.1e}")

# permutation null (bidirectional Bonferroni graph) — EXACT enumeration, matching 07_staging.py.
# The null's support is the set of DISTINCT stage-label assignments to the nodes; it is small
# (1,980 labelings for 12 nodes at composition 1/8/3), so it is enumerated exhaustively rather than
# Monte-Carlo sampled. This makes the plotted null the EXACT null and the printed P the EXACT tail,
# with no seed- or node-order dependence (the earlier 10,000-draw MC estimate wobbled between 1.1e-3
# and 1.2e-3 across node orders for this same edge set). fig1_permnull.csv now holds the exact null
# distribution (one concordance per distinct labeling).
from collections import Counter as _Counter
edges=[(r["exp"],r["out"]) for r in sig]
gnodes=sorted(set(n for e in edges for n in e))
def _msperm(counts):
    n=sum(counts.values()); cur=[None]*n
    def rec(i):
        if i==n:
            yield tuple(cur); return
        for v in sorted(counts):
            if counts[v]:
                counts[v]-=1; cur[i]=v
                yield from rec(i+1)
                counts[v]+=1
    return rec(0)
null=[]
for pm in _msperm(dict(_Counter(STAGE[n] for n in gnodes))):
    smap=dict(zip(gnodes,pm))
    cc=[(a,b) for a,b in edges if smap[a]!=smap[b]]
    ff=sum(1 for a,b in cc if smap[a]<smap[b])
    null.append(ff/len(cc) if cc else 0)
ge=sum(1 for v in null if v>=conc-1e-12); pval=ge/len(null)   # EXACT tail, not add-one
with io.open(os.path.join(OUT,"fig1_permnull.csv"),"w",encoding="utf-8",newline="") as f:
    # full precision, not 5 dp: the published null must let a reader RECOMPUTE the exact tail.
    w=csv.writer(f); w.writerow(["concordance"]); [w.writerow([repr(v)]) for v in null]
with io.open(os.path.join(OUT,"fig1_staging_stat.csv"),"w",encoding="utf-8",newline="") as f:
    w=csv.writer(f); w.writerow(["metric","value"])
    # perm_p is stored at FULL precision. Rounding it to 4 dp here wrote 0.0015 for an exact 3/1980 =
    # 0.001515..., so every consumer that recomputed the exact tail from the null vector disagreed with
    # the canonical file in the 5th decimal. The store must carry the number; the figures do the rounding.
    for k,v in [("observed_concordance",repr(conc)),("perm_p",repr(pval)),("n_cross",len(cross)),
                ("n_forward",len(fwd)),("n_backward",len(back)),("null_mean",f"{sum(null)/len(null):.3f}")]:
        w.writerow([k,v])
print(f"permutation p={pval:.4f} (obs {conc:.3f} vs null mean {sum(null)/len(null):.3f})")

# positive controls with 95% CI
pc_want=[("LDL","CAD"),("BMI","CAD"),("SBP","Stroke")]
with io.open(os.path.join(OUT,"fig1_poscontrols.csv"),"w",encoding="utf-8",newline="") as f:
    w=csv.writer(f); w.writerow(["label","exp","out","b","lo","hi","p","nsnp"])
    for e,o in pc_want:
        r=next((x for x in rows if x["exp"]==e and x["out"]==o),None)
        if r: w.writerow([f"{e}->{o}",e,o,f"{r['b']:.4f}",f"{r['b']-1.96*r['se']:.4f}",f"{r['b']+1.96*r['se']:.4f}",f"{r['p']:.3e}",r["nsnp"]])
print("wrote fig1_edges / fig1_permnull / fig1_staging_stat / fig1_poscontrols")
