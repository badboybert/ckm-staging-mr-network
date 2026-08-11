# -*- coding: utf-8 -*-
"""Step 22: extract TPMI OUTCOME stats at each exposure's clumped instruments.
Exposure & outcome are BOTH TPMI (hg38, same imputation) -> match by MarkerID (chr_pos_A1_A2).
Writes data/harmonised_tpmi/<exp>__<out>.outcome.tsv keyed by the exposure's rsID (so 06-style MR
can harmonise by rsID). Effect allele = Allele2 (as in the exposure). Run: python 22_tpmi_outcomes.py
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

import gzip, io, os, sys
ROOT = CKM_ROOT
BASE = os.path.join(ROOT, "paper 4/independent_build")
UP   = os.path.join(BASE, "data/uploads/tpmi")
INSTR= os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_tpmi"); os.makedirs(HARM, exist_ok=True)

EXPS = ["BMI","SBP","LDL","HDL","TC","TG","HbA1c"]
OUTS = {"T2D":"250.2","CAD":"411.4","HF":"428","Stroke":"433.21","CKD":"585.3"}  # binary PheCodes

def clumped(node):
    f=os.path.join(INSTR,f"{node}.tpmi.clumped.tsv")
    if not os.path.exists(f): return {}
    d={}
    with io.open(f,encoding="utf-8") as fh:
        hdr=fh.readline().rstrip("\n").split("\t"); i_rs=hdr.index("SNP"); i_mid=hdr.index("MarkerID")
        for ln in fh:
            t=ln.rstrip("\n").split("\t"); d[t[i_mid]]=t[i_rs]   # MarkerID -> rsid
    return d

# union of all instrument MarkerIDs (to stream each outcome once)
allmid={}
for e in EXPS:
    for mid,rs in clumped(e).items(): allmid.setdefault(mid, rs)
print(f"[tpmi] union instrument MarkerIDs: {len(allmid)}")

for out,phe in OUTS.items():
    path=os.path.join(UP, f"{phe}.saige.txt.gz"); rows={}
    with gzip.open(path,"rt",errors="replace") as f:
        hdr=f.readline().rstrip("\n").split("\t")
        ix={c:hdr.index(c) for c in ["MarkerID","Allele1","Allele2","AF_Allele2","BETA","SE","p.value","QC"] if c in hdr}
        # binary SAIGE has N_case/N_ctrl; sum for N
        nc = hdr.index("N_case") if "N_case" in hdr else None
        nk = hdr.index("N_ctrl") if "N_ctrl" in hdr else None
        for line in f:
            t=line.rstrip("\n").split("\t")
            try:
                mid=t[ix["MarkerID"]]
                if mid not in allmid: continue
                if t[ix["QC"]]!="PASS": continue
                beta=float(t[ix["BETA"]]); se=float(t[ix["SE"]])
            except (ValueError,IndexError): continue
            if se<=0: continue
            N = (int(float(t[nc]))+int(float(t[nk]))) if (nc is not None and nk is not None) else "NA"
            rows[mid]=(t[ix["Allele2"]].upper(), t[ix["Allele1"]].upper(), t[ix["AF_Allele2"]], beta, se, t[ix["p.value"]], N)
    print(f"[tpmi outcome] {out} ({phe}): matched {len(rows)}/{len(allmid)}")
    for e in EXPS:
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
