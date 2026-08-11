# -*- coding: utf-8 -*-
"""EAS outcome extraction for the hg19-clean subnetwork (BMI/CAD/T2D/eGFR)."""
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

import os, io, sys
sys.path.insert(0, os.path.dirname(__file__))
from lib_formats_eas import EAS_FORMATS, ROOT, _open, load_pos2rs
BASE=os.path.join(ROOT,"paper 4/independent_build"); INSTR=os.path.join(BASE,"data/instruments")
HARM=os.path.join(BASE,"data/harmonised_eas"); os.makedirs(HARM,exist_ok=True)
EXPS=["BMI","SBP","HDL","LDL","TC","TG","HbA1c","CAD","T2D","eGFR"]
OUTS=["CAD","T2D","eGFR","BMI","HDL","LDL","TC","TG","HbA1c"]  # rsID/mappable EAS outcomes (HF/Stroke via 05b)
EDGES=[(e,o) for e in EXPS for o in OUTS if e!=o]

def clumped(node):
    f=os.path.join(INSTR,f"{node}.eas.clumped.tsv")
    if not os.path.exists(f): return set()
    with io.open(f,encoding="utf-8") as fh: fh.readline(); return {ln.split("\t")[0] for ln in fh}

need={}
for e,o in EDGES: need.setdefault(o,set()).update(clumped(e))
pos2rs=None
if any(EAS_FORMATS[o]["mode"]=="map" for o in need): pos2rs=load_pos2rs()

for o,snps in need.items():
    cfg=EAS_FORMATS[o]; C=cfg["cols"]; path=os.path.join(ROOT,cfg["path"]); n_const=cfg.get("n_const")
    if not os.path.exists(path):
        print(f"[EAS outcome] {o}: file not present yet, skipping",flush=True); continue
    rows={}
    with _open(path) as f:
        hdr=f.readline().rstrip("\n").split("\t"); idx={k:hdr.index(v) for k,v in C.items() if v in hdr}
        for line in f:
            t=line.rstrip("\n").split("\t")
            if cfg["mode"]=="rsid":
                try: snp=t[idx["snp"]].strip()
                except IndexError: continue
            else:
                try: snp=pos2rs.get((t[idx["chr"]].strip(),t[idx["pos"]].strip()))
                except IndexError: continue
            if snp not in snps: continue
            try: beta=float(t[idx["beta"]]); se=float(t[idx["se"]])
            except (ValueError,IndexError): continue
            eaf=t[idx["eaf"]] if "eaf" in idx else "NA"
            N=t[idx["n"]] if "n" in idx else (str(n_const) if n_const else "NA")
            rows[snp]=(t[idx["ea"]].strip(),t[idx["oa"]].strip(),eaf,beta,se,t[idx["p"]],N)
    print(f"[EAS outcome] {o}: matched {len(rows)}/{len(snps)}",flush=True)
    for e,oo in EDGES:
        if oo!=o: continue
        with io.open(os.path.join(HARM,f"{e}__{o}.outcome.tsv"),"w",encoding="utf-8") as w:
            w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
            for s in clumped(e):
                if s in rows:
                    ea,oa,eaf,b,se,p,N=rows[s]; w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{o}_EAS\n")
print("done")
