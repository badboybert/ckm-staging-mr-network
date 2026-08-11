# -*- coding: utf-8 -*-
"""Less-conservative TIERED-evidence view of the network (per PI directive): keep nuanced
signals with a graded confidence tier instead of a hard Bonferroni cut. Adds BH-FDR, an
evidence tier, and AHA stage-direction flags. Writes results/network_tiered.csv."""
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

import csv, io, os
BASE = P4_BASE
STAGE = {"BMI":1,"SBP":2,"TG":2,"HDL":2,"TC":2,"LDL":2,"HbA1c":2,"FI":2,"T2D":2,"eGFR":2,"CKD":2,
         "NAFLD":2,"CAD":4,"HF":4,"Stroke":4,"AF":4}

rows=[]
with io.open(os.path.join(BASE,"results/forward_local_edges.csv"),encoding="utf-8") as f:
    for r in csv.DictReader(f):
        try: r["p"]=float(r["ivw_p"]); r["b"]=float(r["ivw_b"])
        except: continue
        r["st"]=str(r["steiger_correct"]).upper()=="TRUE"
        rows.append(r)
rows.sort(key=lambda r:r["p"]); m=len(rows); bonf=0.05/m
for i,r in enumerate(rows,1): r["fdr"]=min(1.0,r["p"]*m/i)
mn=1.0
for r in reversed(rows):
    mn=min(mn,r["fdr"]); r["fdr"]=mn

def tier(r):
    if not r["st"]: return "Steiger-fail"
    if r["p"]<bonf: return "1-Bonferroni"
    if r["fdr"]<0.05: return "2-FDR05"
    if r["fdr"]<0.10: return "3-FDR10"
    if r["p"]<0.05: return "4-nominal"
    return "5-NS"
def sdir(r):
    a,b=STAGE.get(r["exposure"],9),STAGE.get(r["outcome"],9)
    return "forward" if a<b else ("backward" if a>b else "same-stage")

out=os.path.join(BASE,"results/network_tiered.csv")
with io.open(out,"w",encoding="utf-8",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["exposure","outcome","nsnp","ivw_b","ivw_se","ivw_p","fdr","tier","stage_dir","steiger","egger_int_p"])
    for r in rows:
        w.writerow([r["exposure"],r["outcome"],r["nsnp"],f"{r['b']:.4f}",r.get("ivw_se",""),
                    f"{r['p']:.2e}",f"{r['fdr']:.3g}",tier(r),sdir(r),r["steiger_correct"],r.get("egger_intercept_p","")])

# summary
from collections import Counter
tc=Counter(tier(r) for r in rows)
print(f"Tiered network ({m} edges, Bonferroni={bonf:.2e}):")
for t in ["1-Bonferroni","2-FDR05","3-FDR10","4-nominal","5-NS","Steiger-fail"]:
    print(f"  {t:14s}: {tc.get(t,0)}")
print("\nExploratory FEEDBACK/reverse leads (Steiger-correct, backward AHA stage, FDR<0.10):")
for r in rows:
    if r["st"] and sdir(r)=="backward" and r["fdr"]<0.10:
        print(f"  {r['exposure']:6s}->{r['outcome']:6s}  b={r['b']:+.3f} p={r['p']:.2e} fdr={r['fdr']:.3f}  [{tier(r)}]")
print(f"\n-> {out}")
