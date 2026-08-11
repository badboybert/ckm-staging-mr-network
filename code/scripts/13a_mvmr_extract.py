# -*- coding: utf-8 -*-
"""Component A (MVMR) — build harmonised multivariable instrument matrices from LOCAL data.
Union-clump instruments across exposures (PLINK, EUR panel), look up every instrument's
effect on every exposure + outcome in the FULL GWAS, harmonise to the panel reference
allele, drop palindromes/mismatches, write one matrix per model for 13b_mvmr.R.
Run: python 13a_mvmr_extract.py
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

import gzip, io, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_formats import FORMATS, ROOT, _open

BASE  = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
MVDIR = os.path.join(BASE, "data/mvmr"); os.makedirs(MVDIR, exist_ok=True)
PLINK = os.path.join(ROOT, "lhcmr_eas_refs/tools/plink.exe")
PANEL = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur")
BIM   = PANEL + ".bim"

# nodes we need effect-lookups for (exposures + outcomes across all models)
NODES = ["BMI","T2D","SBP","LDL","CAD","HF","Stroke"]
EXPO_FOR_UNION = ["BMI","T2D","SBP","LDL","CAD"]   # instrument sources
MODELS = [
    ("CAD_full",   ["BMI","T2D","LDL","SBP"], "CAD"),
    ("CAD_BMI_T2D",["BMI","T2D"],             "CAD"),
    ("HF_via_CAD", ["BMI","T2D","CAD"],       "HF"),
    ("Stroke_full",["SBP","BMI","LDL"],       "Stroke"),
]

# SBP full GWAS is chr:pos-keyed (Evangelou); handle specially
SBP_RAW = os.path.join(BASE, "data/raw/SBP_Evangelou2018_EUR.txt.gz")
SBP_N   = 757601

def read_clumped(node):
    """rsid -> pval from a clumped instrument file."""
    d={}
    with io.open(os.path.join(INSTR, node+".clumped.tsv"), encoding="utf-8") as f:
        h=f.readline().rstrip("\n").split("\t"); ip=h.index("pval"); isnp=h.index("SNP")
        for line in f:
            t=line.rstrip("\n").split("\t")
            try: d[t[isnp]]=min(d.get(t[isnp],1.0), float(t[ip]))
            except (ValueError,IndexError): continue
    return d

def main():
    # 1. union instruments + min-p
    minp={}
    for node in EXPO_FOR_UNION:
        for rs,p in read_clumped(node).items():
            minp[rs]=min(minp.get(rs,1.0), p)
    print(f"[union] {len(minp)} instruments across {EXPO_FOR_UNION}")

    # 2. PLINK clump the union
    assoc=os.path.join(MVDIR,"union_assoc.txt")
    with io.open(assoc,"w",encoding="utf-8") as w:
        w.write("SNP\tP\n")
        for rs,p in minp.items(): w.write(f"{rs}\t{p}\n")
    out_prefix=os.path.join(MVDIR,"union_clump")
    cmd=[PLINK,"--bfile",PANEL,"--clump",assoc,"--clump-p1","1","--clump-r2","0.001",
         "--clump-kb","10000","--clump-snp-field","SNP","--clump-field","P","--out",out_prefix]
    r=subprocess.run(cmd, capture_output=True, text=True)
    clumped_file=out_prefix+".clumped"
    if not os.path.exists(clumped_file):
        print("[PLINK] clump failed:\n", r.stdout[-1500:], r.stderr[-1500:]); sys.exit(1)
    clumped=set()
    with io.open(clumped_file,encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts=line.split()
            if len(parts)>=3 and parts[2].startswith("rs"): clumped.add(parts[2])
    print(f"[clump] {len(clumped)} LD-independent union instruments")

    # 3. panel bim alleles + positions for the clumped set (reference for harmonisation)
    ref={}   # rsid -> (chr,pos,a1,a2)  ; a1 = reference effect allele
    pos2rs={}
    with io.open(BIM,encoding="utf-8",errors="replace") as f:
        for line in f:
            p=line.split()
            if len(p)<6: continue
            rs=p[1]
            if rs in clumped:
                ref[rs]=(p[0],p[3],p[4].upper(),p[5].upper())
                pos2rs[(p[0],p[3])]=rs
    print(f"[bim] resolved {len(ref)}/{len(clumped)} instruments to panel alleles")

    def palindromic(a1,a2):
        s={a1,a2}; return s=={"A","T"} or s=={"C","G"}

    # 4. extract each node's effect at the clumped instruments
    def extract_standard(node):
        cfg=FORMATS[node]; path=os.path.join(ROOT,cfg["path"]); C=cfg["cols"]
        want=set(ref.keys()); got={}
        with _open(path) as f:
            hdr=f.readline().rstrip("\n").split("\t")
            idx={k:hdr.index(v) for k,v in C.items() if v is not None and v in hdr}
            sep_tab=True
            for line in f:
                t=line.rstrip("\n").split("\t")
                if len(t)<=idx.get("snp",0):
                    t=line.split(); sep_tab=False
                try: rs=t[idx["snp"]].strip()
                except (IndexError,): continue
                if rs not in want: continue
                try:
                    ea=t[idx["ea"]].upper(); oa=t[idx["oa"]].upper()
                    beta=float(t[idx["beta"]]); se=float(t[idx["se"]])
                    eaf=float(t[idx["eaf"]]) if "eaf" in idx else float("nan")
                except (ValueError,IndexError): continue
                if se<=0: continue
                got[rs]=(ea,oa,eaf,beta,se)
                if len(got)==len(want): break
        return got

    def extract_sbp():
        want_pos={(c,p):rs for rs,(c,p,a1,a2) in ref.items()}
        got={}
        with gzip.open(SBP_RAW,"rt",errors="replace") as f:
            hdr=f.readline().split()
            ix={c:hdr.index(c) for c in ["MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","P"]}
            for line in f:
                t=line.split()
                try:
                    mk=t[ix["MarkerName"]].split(":")
                    key=(mk[0],mk[1])
                except (IndexError,): continue
                if key not in want_pos: continue
                rs=want_pos[key]
                try:
                    ea=t[ix["Allele1"]].upper(); oa=t[ix["Allele2"]].upper()
                    beta=float(t[ix["Effect"]]); se=float(t[ix["StdErr"]]); eaf=float(t[ix["Freq1"]])
                except (ValueError,IndexError): continue
                if se<=0: continue
                got[rs]=(ea,oa,eaf,beta,se)
        return got

    cache={}
    for node in NODES:
        cache[node]= extract_sbp() if node=="SBP" else extract_standard(node)
        print(f"[extract] {node}: {len(cache[node])}/{len(ref)} instruments found")

    # 5. harmonise a node's estimate to the panel reference allele a1
    def aligned(node, rs):
        a1,a2=ref[rs][2],ref[rs][3]
        rec=cache[node].get(rs)
        if rec is None: return None
        ea,oa,eaf,beta,se=rec
        if {ea,oa}!={a1,a2}: return None          # allele mismatch (indel/strand) -> drop
        if ea==a1: return (beta,se,eaf)
        else:      return (-beta,se,(1-eaf) if eaf==eaf else eaf)  # flip to a1

    # 6. assemble per-model matrices
    for name,expos,outc in MODELS:
        need=expos+[outc]
        rows=[]
        for rs,(c,p,a1,a2) in ref.items():
            if palindromic(a1,a2): continue
            vals={}
            ok=True
            for node in need:
                a=aligned(node,rs)
                if a is None: ok=False; break
                vals[node]=a
            if not ok: continue
            rows.append((rs,vals))
        # write
        out=os.path.join(MVDIR,f"model_{name}.tsv")
        with io.open(out,"w",encoding="utf-8",newline="") as w:
            cols=["SNP"]
            for e in expos: cols += [f"bx_{e}", f"se_{e}"]
            cols += ["by","sey"]
            w.write("\t".join(cols)+"\n")
            for rs,vals in rows:
                line=[rs]
                for e in expos: line += [f"{vals[e][0]:.6g}", f"{vals[e][1]:.6g}"]
                line += [f"{vals[outc][0]:.6g}", f"{vals[outc][1]:.6g}"]
                w.write("\t".join(line)+"\n")
        print(f"[model] {name}: {outc} ~ {'+'.join(expos)}  -> {len(rows)} SNPs  ({os.path.basename(out)})")

if __name__=="__main__":
    main()
