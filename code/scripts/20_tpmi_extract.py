# -*- coding: utf-8 -*-
"""TPMI (Taiwan, EAS) SAIGE ingest — extract genome-wide-sig instruments.
TPMI is hg38, chr:pos-keyed (MarkerID=chr_pos_Allele1_Allele2; effect allele = Allele2, AF_Allele2/BETA
per Allele2), with a QC flag + imputationInfo. Route to rsID (for clumping): TPMI hg38 pos -> liftover
hg38->hg19 -> rsID via the 1000G-EAS pos2rs cache. Exposure<->outcome later match by MarkerID (both TPMI).
Run: python 20_tpmi_extract.py [nodes...]
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

import gzip, io, os, sys, pickle
from pyliftover import LiftOver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = CKM_ROOT
BASE = os.path.join(ROOT, "paper 4/independent_build")
UP   = os.path.join(BASE, "data/uploads/tpmi")
INSTR= os.path.join(BASE, "data/instruments"); os.makedirs(INSTR, exist_ok=True)
CACHE= os.path.join(BASE, "data/eas_pos2rs.pkl")   # hg19 (chr,pos)->rsid from merged EAS bim

# node -> filename in uploads/tpmi/ . SAIGE cols are fixed; effect allele = Allele2.
TPMI = {
 # exposures (quantitative)
 "BMI":"BMI","SBP":"SBP","DBP":"DBP","LDL":"LDL-C","HDL":"HDL-C","TC":"TC","TG":"TG",
 "HbA1c":"HbA1c","FG":"FG","sCr":"sCR","WBC":"WBC","ALT":"ALT","AST":"AST",
 # outcomes (binary, PheCodes)
 "T2D":"250.2","CAD":"411.4","HF":"428","Stroke":"433.21","AF":"427.21","MASLD":"571.5","CKD":"585.3",
}

def fpath(node): return os.path.join(UP, f"{TPMI[node]}.saige.txt.gz")

def load_pos2rs(log=print):
    log("[tpmi] loading hg19 pos2rs cache...")
    with open(CACHE,"rb") as f: d=pickle.load(f)
    log(f"[tpmi] pos2rs entries: {len(d):,}")
    return d

def extract_sig(node, pthr=5e-8, info_min=0.3, log=print):
    """stream TPMI SAIGE, keep QC=PASS & info>=info_min & p<pthr; return list of sig rows (hg38)."""
    path=fpath(node); rows=[]; seen=0
    with gzip.open(path,"rt",errors="replace") as f:
        hdr=f.readline().rstrip("\n").split("\t")
        # QT files have 'N'; binary (PheCode) files have 'N_case'+'N_ctrl' instead
        need=["CHR","POS","MarkerID","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","QC"]
        ix={c:hdr.index(c) for c in need if c in hdr}
        i_N = hdr.index("N") if "N" in hdr else None
        i_nc = hdr.index("N_case") if "N_case" in hdr else None
        i_nk = hdr.index("N_ctrl") if "N_ctrl" in hdr else None
        for line in f:
            seen+=1; t=line.rstrip("\n").split("\t")
            try:
                if t[ix["QC"]]!="PASS": continue
                if float(t[ix["imputationInfo"]])<info_min: continue
                p=float(t[ix["p.value"]])
            except (ValueError,IndexError): continue
            if not (p<pthr): continue
            try:
                chrom=t[ix["CHR"]].replace("chr",""); pos=int(t[ix["POS"]])
                ea=t[ix["Allele2"]].upper(); oa=t[ix["Allele1"]].upper()
                eaf=float(t[ix["AF_Allele2"]]); beta=float(t[ix["BETA"]]); se=float(t[ix["SE"]])
                mid=t[ix["MarkerID"]]
                N = t[i_N] if i_N is not None else (str(int(float(t[i_nc]))+int(float(t[i_nk]))) if (i_nc is not None and i_nk is not None) else "NA")
            except (ValueError,IndexError): continue
            if se<=0: continue
            # skip multiallelic-ambiguous long indels for clumping robustness (keep SNVs + short indels)
            rows.append((mid,chrom,pos,ea,oa,eaf,beta,se,p,N))
    log(f"[tpmi] {node}: {len(rows)} sig (p<{pthr:g}, PASS, info>={info_min}) from {seen:,} rows")
    return rows

def main():
    nodes = sys.argv[1:] or ["BMI","SBP","LDL","HDL","TC","TG","HbA1c"]
    lo = LiftOver('hg38','hg19')
    pos2rs = load_pos2rs()
    for node in nodes:
        rows = extract_sig(node)
        out = os.path.join(INSTR, f"{node}.tpmi.sig.tsv")
        kept=0; nomap=0
        with io.open(out,"w",encoding="utf-8") as w:
            w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tMarkerID\tPhenotype\n")
            for mid,chrom,pos,ea,oa,eaf,beta,se,p,N in rows:
                conv = lo.convert_coordinate(f"chr{chrom}", pos-1)   # pyliftover is 0-based
                if not conv: nomap+=1; continue
                hg19pos = str(conv[0][1]+1)   # pos2rs keys store position as a STRING
                rs = pos2rs.get((chrom, hg19pos)) or pos2rs.get((f"chr{chrom}", hg19pos))
                if not rs or not str(rs).startswith("rs"): nomap+=1; continue
                w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{beta}\t{se}\t{p}\t{N}\t{mid}\t{node}\n"); kept+=1
        print(f"[tpmi] {node}: {kept} mapped to rsID (hg38->hg19->rs), {nomap} unmapped -> {os.path.basename(out)}", flush=True)

if __name__=="__main__": main()
