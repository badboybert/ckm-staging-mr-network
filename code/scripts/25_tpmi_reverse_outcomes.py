# -*- coding: utf-8 -*-
"""Step 25: TPMI REVERSE edges — extract RF-as-OUTCOME stats at each disease's clumped instruments.
Mirror of step 22 with exposures=diseases, outcomes=risk-factors. Match by MarkerID (both TPMI).
Writes data/harmonised_tpmi/<disease>__<RF>.outcome.tsv. Run: python 25_tpmi_reverse_outcomes.py
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

import gzip, io, os
ROOT = CKM_ROOT
BASE = os.path.join(ROOT, "paper 4/independent_build")
UP   = os.path.join(BASE, "data/uploads/tpmi")
INSTR= os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_tpmi"); os.makedirs(HARM, exist_ok=True)

DIS  = ["T2D","CAD","HF","Stroke","CKD"]                 # exposures (binary)
RF   = {"BMI":"BMI","SBP":"SBP","LDL":"LDL-C","HDL":"HDL-C","TC":"TC","TG":"TG","HbA1c":"HbA1c"}  # outcomes (quant)

def clumped(node):
    f=os.path.join(INSTR,f"{node}.tpmi.clumped.tsv")
    if not os.path.exists(f): return {}
    d={}
    with io.open(f,encoding="utf-8") as fh:
        hdr=fh.readline().rstrip("\n").split("\t"); i_rs=hdr.index("SNP"); i_mid=hdr.index("MarkerID")
        for ln in fh:
            t=ln.rstrip("\n").split("\t"); d[t[i_mid]]=t[i_rs]
    return d

allmid={}
for e in DIS:
    for mid,rs in clumped(e).items(): allmid.setdefault(mid, rs)
print(f"[tpmi rev] union disease-instrument MarkerIDs: {len(allmid)}")

for out,phe in RF.items():
    path=os.path.join(UP, f"{phe}.saige.txt.gz"); rows={}
    with gzip.open(path,"rt",errors="replace") as f:
        hdr=f.readline().rstrip("\n").split("\t")
        ix={c:hdr.index(c) for c in ["MarkerID","Allele1","Allele2","AF_Allele2","BETA","SE","p.value","N","QC"] if c in hdr}
        for line in f:
            t=line.rstrip("\n").split("\t")
            try:
                mid=t[ix["MarkerID"]]
                if mid not in allmid: continue
                if t[ix["QC"]]!="PASS": continue
                beta=float(t[ix["BETA"]]); se=float(t[ix["SE"]])
            except (ValueError,IndexError): continue
            if se<=0: continue
            N = t[ix["N"]] if "N" in ix else "NA"
            rows[mid]=(t[ix["Allele2"]].upper(), t[ix["Allele1"]].upper(), t[ix["AF_Allele2"]], beta, se, t[ix["p.value"]], N)
    print(f"[tpmi rev outcome] {out} ({phe}): matched {len(rows)}/{len(allmid)}")
    for e in DIS:
        exp_mid=clumped(e)
        outf=os.path.join(HARM, f"{e}__{out}.outcome.tsv"); m=0
        with io.open(outf,"w",encoding="utf-8") as w:
            w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
            for mid,rs in exp_mid.items():
                if mid in rows:
                    ea,oa,eaf,b,se,p,N=rows[mid]
                    w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{out}\n"); m+=1
        print(f"   {e}->{out}: {m} outcome SNPs")
print("done")
