# -*- coding: utf-8 -*-
"""Extract hg38 EAS outcomes (HF/Stroke GCST90-series) for hg19 instrument rsIDs.
Route: exposure rsID -> hg19 (chr,pos) from EAS bim -> liftover hg19->hg38 -> match hg38 outcome by pos.
Only the ~hundreds of instrument SNPs are lifted (cheap)."""
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

import os, io, gzip, sys
from pyliftover import LiftOver
sys.path.insert(0, os.path.dirname(__file__))
from lib_formats_eas import ROOT, BIM

BASE=os.path.join(ROOT,"paper 4/independent_build"); INSTR=os.path.join(BASE,"data/instruments")
HARM=os.path.join(BASE,"data/harmonised_eas"); os.makedirs(HARM,exist_ok=True)

# hg38 EAS outcomes: file, plain/gz, col map (chr,pos,ea,oa,eaf,beta,se,p[,n])
HG38 = {
 "HF": dict(path="paper 6/analysis/data/outcomes/eas/HF_GCST90668009.gz", gz=False,
            cols=dict(chr="chromosome",pos="base_pair_location",ea="effect_allele",oa="other_allele",
                      eaf="effect_allele_frequency",beta="beta",se="standard_error",p="p_value",n="n")),
 "Stroke": dict(path="paper 6/analysis/data/outcomes/eas/Stroke_EAS_GCST90104545.tsv.gz", gz=True,
            cols=dict(chr="chromosome",pos="base_pair_location",ea="effect_allele",oa="other_allele",
                      eaf="effect_allele_frequency",beta="beta",se="standard_error",p="p_value"), n_const=256274),  # GIGASTROKE EAS ischemic stroke (19,032 cases + 237,242 controls); was mislabeled 1,296,908
}
EXPS=["BMI","SBP","HDL","LDL","TC","TG","HbA1c","CAD","T2D","eGFR"]

def clumped(node):
    f=os.path.join(INSTR,f"{node}.eas.clumped.tsv")
    if not os.path.exists(f): return []
    with io.open(f,encoding="utf-8") as fh: fh.readline(); return [ln.split("\t")[0] for ln in fh]

# 1. union of all instrument rsIDs
allrs=set()
for e in EXPS: allrs.update(clumped(e))
print(f"instrument rsIDs to place: {len(allrs)}",flush=True)

# 2. rsID -> hg19 (chr,pos) from EAS bim
rs2hg19={}
with io.open(BIM,errors="replace") as f:
    for ln in f:
        c,rs,cm,bp,a1,a2=ln.split()
        if rs in allrs: rs2hg19[rs]=(c,int(bp))
print(f"placed on hg19: {len(rs2hg19)}",flush=True)

# 3. build a COMBINED position->rsid lookup covering BOTH hg19 and hg38 positions,
#    so we match regardless of the outcome file's build (GCST90 EAS files vary).
lo=LiftOver('hg19','hg38')
pos2rs={}
for rs,(c,bp) in rs2hg19.items():
    pos2rs[(c, bp)]=rs                       # hg19
    r=lo.convert_coordinate('chr'+c, bp)
    if r: pos2rs[(c, r[0][1])]=rs            # hg38 (if outcome is hg38)
hg38pos2rs=pos2rs
print(f"combined hg19+hg38 lookup keys: {len(pos2rs)}",flush=True)

# 4. stream each hg38 outcome, match by (chr,hg38pos)
def _open(p,gz): return gzip.open(p,"rt",errors="replace") if gz else open(p,"rt",errors="replace")
for onode,cfg in HG38.items():
    C=cfg["cols"]; rows={}; n_const=cfg.get("n_const")
    with _open(os.path.join(ROOT,cfg["path"]), cfg["gz"]) as f:
        hdr=f.readline().rstrip("\n").split("\t"); idx={k:hdr.index(v) for k,v in C.items() if v in hdr}
        for line in f:
            t=line.rstrip("\n").split("\t")
            try: key=(t[idx["chr"]].strip(), int(t[idx["pos"]]))
            except (ValueError,IndexError): continue
            rs=hg38pos2rs.get(key)
            if not rs: continue
            try: beta=float(t[idx["beta"]]); se=float(t[idx["se"]])
            except (ValueError,IndexError): continue
            eaf=t[idx["eaf"]] if "eaf" in idx else "NA"
            N=t[idx["n"]] if "n" in idx else (str(n_const) if n_const else "NA")
            rows[rs]=(t[idx["ea"]].strip(),t[idx["oa"]].strip(),eaf,beta,se,t[idx["p"]],N)
    print(f"[hg38 outcome] {onode}: matched {len(rows)}",flush=True)
    for e in EXPS:
        with io.open(os.path.join(HARM,f"{e}__{onode}.outcome.tsv"),"w",encoding="utf-8") as w:
            w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
            for s in clumped(e):
                if s in rows:
                    ea,oa,eaf,b,se,p,N=rows[s]; w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{onode}_EAS\n")
print("done")
