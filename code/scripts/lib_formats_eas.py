# -*- coding: utf-8 -*-
"""EAS format registry + chr:pos->rsid mapper (via merged hg19 1000G-EAS bim).
Handles the two EAS file shapes: rsID-present (clump directly) and chr:pos-only (map).
BUILD NOTE: BMI(Akiyama)/CAD(BBJ)/T2D(AGEN)/eGFR(BBJ hum0014) are hg19 = matches the panel.
The GCST90-series EAS outcomes (HF/Stroke/CKD/ALT) are hg38 -> need liftover before chr:pos map (deferred)."""
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

import gzip, io, os, math, pickle
ROOT = CKM_ROOT
BIM  = os.path.join(ROOT, "paper 4/independent_build/data/g1000_eas_merged.bim")
CACHE= os.path.join(ROOT, "paper 4/independent_build/data/eas_pos2rs.pkl")

# node -> path, build, snp-mode ('rsid' col or 'map' chr:pos), column map
EAS_FORMATS = {
 "BMI":  dict(path="paper 4/independent_build/data/raw/BMI_BBJ_Akiyama2017_EAS.txt.gz", build="hg19", mode="map",
              cols=dict(chr="CHR",pos="POS",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P"), n_const=158284),
 "CAD":  dict(path="paper 3.eas/analysis/data/eas_outcomes/BBJ_CAD.txt.gz", build="hg19", mode="rsid",
              cols=dict(snp="SNPID",ea="Allele2",oa="Allele1",eaf="AF_Allele2",beta="BETA",se="SE",p="p.value",n="N")),
 "T2D":  dict(path="paper 3.eas/analysis/data/eas_outcomes/AGEN_T2D_full.txt.gz", build="hg19", mode="map",
              cols=dict(chr="Chr",pos="Pos",ea="EA",oa="NEA",eaf="EAF",beta="Beta",se="SE",p="P",n="Neff")),
 "eGFR": dict(path="paper 3.eas/analysis/data/eas_outcomes/hum0014.v7.eGFR/BBJ.eGFR.autosome.txt", build="hg19", mode="rsid",
              cols=dict(snp="SNP",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P",n="N")),
 # --- BBJ Kanai2018 QT exposures from NBDC hum0014 (rsID, hg19, same BBJ.<trait>.autosome format) ---
 "SBP":  dict(path="paper 4/independent_build/data/raw/hum0014.v7.SBP/BBJ.SBP.autosome.txt", build="hg19", mode="rsid",
              cols=dict(snp="SNP",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P",n="N")),
 "HbA1c":dict(path="paper 4/independent_build/data/raw/hum0014.v7.HbA1c/BBJ.HbA1c.autosome.txt", build="hg19", mode="rsid",
              cols=dict(snp="SNP",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P",n="N")),
 "HDL":  dict(path="paper 4/independent_build/data/raw/hum0014.v7.HDL/BBJ.HDL-C.autosome.txt", build="hg19", mode="rsid",
              cols=dict(snp="SNP",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P",n="N")),
 "LDL":  dict(path="paper 4/independent_build/data/raw/hum0014.v7.LDL/BBJ.LDL-C.autosome.txt", build="hg19", mode="rsid",
              cols=dict(snp="SNP",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P",n="N")),
 "TC":   dict(path="paper 4/independent_build/data/raw/hum0014.v7.TC/BBJ.TC.autosome.txt", build="hg19", mode="rsid",
              cols=dict(snp="SNP",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P",n="N")),
 "TG":   dict(path="paper 4/independent_build/data/raw/hum0014.v7.TG/BBJ.TG.autosome.txt", build="hg19", mode="rsid",
              cols=dict(snp="SNP",ea="ALT",oa="REF",eaf="Frq",beta="BETA",se="SE",p="P",n="N")),
 # hg38 GCST90 EAS outcomes — DEFERRED (need hg38->hg19 liftover before chr:pos mapping):
 # HF GCST90668009, Stroke GCST90104545, CKD GCST90018602, ALT GCST90278652
}

def _open(p): return gzip.open(p,"rt",errors="replace") if p.endswith(".gz") else open(p,"rt",errors="replace")

def load_pos2rs(log=print):
    if os.path.exists(CACHE):
        with open(CACHE,"rb") as f: return pickle.load(f)
    log("[EAS] building (chr,pos)->rsid map from merged bim (one-time)...")
    m={}
    with io.open(BIM,errors="replace") as f:
        for ln in f:
            c,rs,cm,bp,a1,a2=ln.split()
            m[(c,bp)]=rs
    with open(CACHE,"wb") as f: pickle.dump(m,f)
    log(f"[EAS] pos2rs map: {len(m)} entries -> cached")
    return m

def extract_sig_eas(node, pthr=5e-8, log=print, pos2rs=None):
    cfg=EAS_FORMATS[node]; C=cfg["cols"]; path=os.path.join(ROOT,cfg["path"])
    out_dir=os.path.join(ROOT,"paper 4/independent_build/data/instruments")
    os.makedirs(out_dir,exist_ok=True); out=os.path.join(out_dir,f"{node}.eas.sig.tsv")
    n_const=cfg.get("n_const")
    if cfg["mode"]=="map" and pos2rs is None: pos2rs=load_pos2rs(log)
    kept=0
    with _open(path) as f, io.open(out,"w",encoding="utf-8") as w:
        hdr=f.readline().rstrip("\n").split("\t")
        idx={k:hdr.index(v) for k,v in C.items() if v in hdr}
        w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
        for line in f:
            t=line.rstrip("\n").split("\t")
            try: p=float(t[idx["p"]])
            except (ValueError,IndexError): continue
            if not (p<pthr): continue
            try:
                ea=t[idx["ea"]].strip(); oa=t[idx["oa"]].strip()
                beta=float(t[idx["beta"]]); se=float(t[idx["se"]])
            except (ValueError,IndexError): continue
            if se<=0 or math.isnan(beta) or math.isnan(se): continue
            if cfg["mode"]=="rsid":
                snp=t[idx["snp"]].strip()
            else:
                snp=pos2rs.get((t[idx["chr"]].strip(),t[idx["pos"]].strip()))
            if not snp or not snp.startswith("rs"): continue
            eaf=t[idx["eaf"]] if "eaf" in idx else "NA"
            N=t[idx["n"]] if "n" in idx else (str(n_const) if n_const else "NA")
            w.write(f"{snp}\t{ea}\t{oa}\t{eaf}\t{beta}\t{se}\t{p}\t{N}\t{node}_EAS\n"); kept+=1
    log(f"[EAS extract] {node}: {kept} sig SNPs -> {os.path.basename(out)}")
    return out

if __name__=="__main__":
    import sys
    p2r=None
    nodes=sys.argv[1:] or ["BMI","CAD","T2D","eGFR"]
    if any(EAS_FORMATS[n]["mode"]=="map" for n in nodes): p2r=load_pos2rs()
    for n in nodes: extract_sig_eas(n, pos2rs=p2r)
