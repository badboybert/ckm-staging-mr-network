# -*- coding: utf-8 -*-
"""Per-file column-format registry + streaming genome-wide-significant instrument extractor.
Standardises heterogeneous on-disk GWAS to: SNP EA OA EAF BETA SE P N (+Phenotype)."""
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

import gzip, io, os, math

ROOT = CKM_ROOT

# node -> dict(path, cols{snp,ea,oa,eaf,beta,se,p,n}, [n_const], [sep])
FORMATS = {
 "CAD":  dict(path="paper 6/analysis/data/outcomes/CAD_GCST90132314.h.tsv.gz",
              cols=dict(snp="rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n="n")),
 "HF":   dict(path="paper 6/analysis/data/outcomes/HF_GCST90728695.tsv.gz",
              cols=dict(snp="rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n="N_total")),
 "T2D":  dict(path="paper 6/analysis/data/outcomes/T2D_GCST006867.h.tsv.gz",
              cols=dict(snp="hm_rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n="n")),
 "LDL":  dict(path="paper 6/analysis/data/outcomes/LDL_GLGC2021_EUR.gz",
              cols=dict(snp="rsID",ea="ALT",oa="REF",eaf="POOLED_ALT_AF",
                        beta="EFFECT_SIZE",se="SE",p="pvalue",n="N")),
 "eGFR": dict(path="paper 6/analysis/data/outcomes/eur_upgraded/eGFR_Stanzick2021_EUR.tsv",
              cols=dict(snp="variant_id",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n="n")),
 "Stroke":dict(path="paper 6/analysis/data/outcomes/Stroke_GCST90104539.h.tsv.gz",
              cols=dict(snp="rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n=None), n_const=1308460),
 # --- HF subtypes (Henry 2024 HERMES; same batch/format as HF above; N_total = cases+controls) ---
 "niHF":   dict(path="paper 4/independent_build/data/raw/HF_niHF_GCST90728696_EUR.tsv.gz",
              cols=dict(snp="rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n="N_total")),
 "niHFpEF":dict(path="paper 4/independent_build/data/raw/HF_niHFpEF_GCST90728697_EUR.tsv.gz",
              cols=dict(snp="rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n="N_total")),
 "niHFrEF":dict(path="paper 4/independent_build/data/raw/HF_niHFrEF_GCST90728698_EUR.tsv.gz",
              cols=dict(snp="rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n="N_total")),
 # --- downloaded EUR metabolic exposures (Step 2) ---
 "BMI":  dict(path="paper 4/independent_build/data/raw/BMI_Yengo2018_EUR.txt.gz",
              cols=dict(snp="SNP",ea="Tested_Allele",oa="Other_Allele",eaf="Freq_Tested_Allele_in_HRS",
                        beta="BETA",se="SE",p="P",n="N")),
 "HbA1c":dict(path="paper 4/independent_build/data/raw/HbA1c_MAGIC_EUR.h.tsv.gz",
              cols=dict(snp="hm_rsid",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n=None), n_const=146806),
 "FI":   dict(path="paper 4/independent_build/data/raw/FI_MAGIC_EUR.f.tsv.gz",
              cols=dict(snp="variant_id",ea="effect_allele",oa="other_allele",eaf="effect_allele_frequency",
                        beta="beta",se="standard_error",p="p_value",n=None), n_const=151013),
 "HDL":  dict(path="paper 4/independent_build/data/raw/HDL_GLGC2021_EUR.gz",
              cols=dict(snp="rsID",ea="ALT",oa="REF",eaf="POOLED_ALT_AF",
                        beta="EFFECT_SIZE",se="SE",p="pvalue",n="N")),
 "TC":   dict(path="paper 4/independent_build/data/raw/TC_GLGC2021_EUR.gz",
              cols=dict(snp="rsID",ea="ALT",oa="REF",eaf="POOLED_ALT_AF",
                        beta="EFFECT_SIZE",se="SE",p="pvalue",n="N")),
 "TG":   dict(path="paper 4/independent_build/data/raw/TG_GLGC2021_EUR.gz",
              cols=dict(snp="rsID",ea="ALT",oa="REF",eaf="POOLED_ALT_AF",
                        beta="EFFECT_SIZE",se="SE",p="pvalue",n="N")),
 # --- ApoB (Sinnott-Armstrong 2021 UKB, GWAS-SSF harmonised; item 1 atherogenic-axis MVMR) ---
 "ApoB": dict(path="paper 4/independent_build/data/uploads/apob/ApoB_SinnottArmstrong2021_GCST90025952_EUR.h.tsv.gz",
              cols=dict(snp="hm_rsid",ea="hm_effect_allele",oa="hm_other_allele",eaf="hm_effect_allele_frequency",
                        beta="hm_beta",se="standard_error",p="p_value",n=None), n_const=435744),
}

def _open(p):
    return gzip.open(p,"rt",errors="replace") if p.endswith(".gz") else open(p,"rt",errors="replace")

def extract_sig(node, pthr=5e-8, out_dir=None, log=print):
    """Stream the node's GWAS, emit standardized rows with p<pthr. Returns output path."""
    cfg = FORMATS[node]; path=os.path.join(ROOT,cfg["path"]); C=cfg["cols"]
    out_dir = out_dir or os.path.join(ROOT,"paper 4/independent_build/data/instruments")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{node}.sig.tsv")
    n_const = cfg.get("n_const")
    kept=0; seen=0
    with _open(path) as f, io.open(out,"w",encoding="utf-8") as w:
        hdr=f.readline().rstrip("\n").split("\t")
        idx={name:hdr.index(col) for name,col in C.items() if col is not None and col in hdr}
        need=["snp","ea","oa","beta","se","p"]
        missing=[k for k in need if k not in idx]
        if missing: raise KeyError(f"{node}: missing columns {missing}; header={hdr[:8]}")
        w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
        for line in f:
            seen+=1
            t=line.rstrip("\n").split("\t")
            try:
                p=float(t[idx["p"]])
            except (ValueError,IndexError):
                continue
            if not (p<pthr): continue
            try:
                snp=t[idx["snp"]].strip(); ea=t[idx["ea"]].strip(); oa=t[idx["oa"]].strip()
                beta=float(t[idx["beta"]]); se=float(t[idx["se"]])
            except (ValueError,IndexError):
                continue
            if not snp.startswith("rs"): continue
            if se<=0 or math.isnan(beta) or math.isnan(se): continue
            eaf = t[idx["eaf"]] if "eaf" in idx else "NA"
            N = (t[idx["n"]] if "n" in idx else str(n_const)) if ("n" in idx or n_const) else "NA"
            w.write(f"{snp}\t{ea}\t{oa}\t{eaf}\t{beta}\t{se}\t{p}\t{N}\t{node}\n")
            kept+=1
    log(f"[extract] {node}: {kept} genome-wide-sig SNPs (p<{pthr:g}) from {seen} rows -> {os.path.basename(out)}")
    return out

if __name__=="__main__":
    import sys
    for node in (sys.argv[1:] or ["LDL","CAD","HF"]):
        extract_sig(node)
