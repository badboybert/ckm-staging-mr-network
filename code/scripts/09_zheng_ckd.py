# -*- coding: utf-8 -*-
"""Zheng2021 replication: SBP->CKD and BMI->CKD in EUR and EAS (binary CKD outcome).
Self-contained: load clumped instruments -> stream CKD outcome (EUR=OR->lnOR by rsid;
EAS=beta, rsid-or-chr:pos combined lookup) -> harmonise alleles -> Python IVW + Egger-intercept."""
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

import io, os, gzip, math, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy import stats

ROOT=CKM_ROOT; BASE=os.path.join(ROOT,"paper 4/independent_build")
INSTR=os.path.join(BASE,"data/instruments")

def load_instr(f):
    d={}
    with io.open(f,encoding="utf-8") as fh:
        fh.readline()
        for ln in fh:
            t=ln.rstrip("\n").split("\t")
            try: d[t[0]]=(t[1].upper(),t[2].upper(),float(t[4]),float(t[5]))  # ea,oa,beta,se
            except: pass
    return d

def ivw(instr, out):
    """out: rs-> (ea,oa,beta,se). Harmonise to exposure effect allele; IVW + Egger intercept."""
    bx=[];by=[];sy=[]
    for rs,(eea,eoa,eb,ese) in instr.items():
        if rs not in out: continue
        oea,ooa,ob,ose=out[rs]
        if {eea,eoa}!={oea,ooa}: continue         # allele mismatch
        if oea==eea: pass
        elif oea==eoa: ob=-ob                       # flip to exposure effect allele
        else: continue
        bx.append(eb); by.append(ob); sy.append(ose)
    if len(bx)<3: return None
    bx=np.array(bx);by=np.array(by);sy=np.array(sy); w=1/sy**2
    b=np.sum(w*bx*by)/np.sum(w*bx**2); se=math.sqrt(1/np.sum(w*bx**2))
    p=2*stats.norm.sf(abs(b/se))
    return dict(n=len(bx), b=b, se=se, p=p)

# CKD outcome loaders
def eur_ckd_lookup(rsset):
    # CKDGen Wuttke2019 EUR CKD (GCST008065): whitespace-delimited; Effect=beta(logOR), StdErr, RSID, Allele1=EA
    p=os.path.join(ROOT,"paper 6/analysis/data/outcomes/CKD_GCST008065.txt.gz"); out={}
    with gzip.open(p,"rt",errors="replace") as f:
        h=f.readline().split(); ix={c:h.index(c) for c in ["RSID","Allele1","Allele2","Effect","StdErr"]}
        for ln in f:
            t=ln.split()
            try:
                rs=t[ix["RSID"]]
                if rs not in rsset: continue
                b=float(t[ix["Effect"]]); se=float(t[ix["StdErr"]])
            except: continue
            out[rs]=(t[ix["Allele1"]].upper(),t[ix["Allele2"]].upper(),b,se)
    return out

def eas_ckd_lookup(rsset, pos2rs_inv):
    """EAS CKD: beta present; variant_id sparse -> match by rsid if present else (chr,pos) via inv map."""
    p=os.path.join(ROOT,"paper 6/analysis/data/outcomes/eas/CKD_GCST90018602.gz"); out={}
    with gzip.open(p,"rt",errors="replace") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={c:h.index(c) for c in ["chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","variant_id"]}
        for ln in f:
            t=ln.rstrip("\n").split("\t")
            rs=t[ix["variant_id"]].strip()
            if not rs.startswith("rs"):
                rs=pos2rs_inv.get((t[ix["chromosome"]].strip(),t[ix["base_pair_location"]].strip()))
            if not rs or rs not in rsset: continue
            try: b=float(t[ix["beta"]]); se=float(t[ix["standard_error"]])
            except: continue
            out[rs]=(t[ix["effect_allele"]].upper(),t[ix["other_allele"]].upper(),b,se)
    return out

# --- EUR: SBP->CKD, BMI->CKD ---
print("=== Zheng2021 replication: SBP->CKD and BMI->CKD ===\n")
for exp in ["SBP","BMI"]:
    ins=load_instr(os.path.join(INSTR,f"{exp}.clumped.tsv"))
    out=eur_ckd_lookup(set(ins))
    r=ivw(ins,out)
    if r: print(f"EUR  {exp}->CKD : nSNP={r['n']:3d}  IVW b={r['b']:+.4f} se={r['se']:.4f} p={r['p']:.2e}")
    else: print(f"EUR  {exp}->CKD : <3 SNPs")

# --- EAS: need pos2rs for chr:pos fallback (build the merged EAS bim inverse) ---
from lib_formats_eas import BIM
pos2rs={}
with io.open(BIM,errors="replace") as f:
    for ln in f:
        c,rs,cm,bp,a1,a2=ln.split(); pos2rs[(c,bp)]=rs
for exp in ["SBP","BMI"]:
    ins=load_instr(os.path.join(INSTR,f"{exp}.eas.clumped.tsv"))
    out=eas_ckd_lookup(set(ins), pos2rs)
    r=ivw(ins,out)
    if r: print(f"EAS  {exp}->CKD : nSNP={r['n']:3d}  IVW b={r['b']:+.4f} se={r['se']:.4f} p={r['p']:.2e}")
    else: print(f"EAS  {exp}->CKD : <3 SNPs (matched {len(out)})")
print("\nNOTE: CKD is binary (log-OR). Zheng2021: HTN/SBP->CKD causal in EUR, NULL in EAS; BMI->CKD both (nonlinear threshold).")
