# -*- coding: utf-8 -*-
"""Close the 5 ledger NO_REVERSE_DATA gaps. Build outcome files for the reverse edges:
  Stroke->{BMI,T2D,CAD,HF}  (Stroke as EXPOSURE; rsID outcomes via lib_formats FORMATS)
  {CAD,HF,Stroke}->SBP      (SBP as OUTCOME; SBP is chr:pos-keyed Evangelou -> rsid->chr:pos via EUR bim)
Writes data/harmonised/<exp>__<out>.outcome.tsv (same schema as 05_extract_outcomes.py).
Run: python 16_reverse_outcomes.py   (after Stroke.clumped.tsv exists)
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

import io, os, sys, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_formats import FORMATS, ROOT, _open

BASE  = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments"); HARM = os.path.join(BASE, "data/harmonised")
BIM   = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur.bim")
SBP_RAW = os.path.join(BASE, "data/raw/SBP_Evangelou2018_EUR.txt.gz"); SBP_N = 757601

def clumped_snps(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f): return []
    with io.open(f, encoding="utf-8") as fh:
        fh.readline(); return [ln.split("\t")[0] for ln in fh if ln.strip()]

def write_outcome(exp, out, rows_by_snp, out_label):
    exp_snps = clumped_snps(exp)
    outf = os.path.join(HARM, f"{exp}__{out}.outcome.tsv")
    m = 0
    with io.open(outf, "w", encoding="utf-8") as w:
        w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
        for s in exp_snps:
            if s in rows_by_snp:
                ea, oa, eaf, b, se, p, N = rows_by_snp[s]
                w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{out_label}\n"); m += 1
    print(f"  {exp}->{out}: {m}/{len(exp_snps)} outcome SNPs -> {os.path.basename(outf)}", flush=True)

# ---- 1. Stroke as EXPOSURE -> rsID outcomes (BMI, T2D, CAD, HF) ----
STROKE_OUTS = ["BMI","T2D","CAD","HF"]
stroke_snps = set(clumped_snps("Stroke"))
print(f"[Stroke exposure] {len(stroke_snps)} clumped instruments")
for out in STROKE_OUTS:
    cfg = FORMATS[out]; C = cfg["cols"]; path = os.path.join(ROOT, cfg["path"]); n_const = cfg.get("n_const")
    rows = {}
    with _open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: hdr.index(v) for k, v in C.items() if v and v in hdr}
        for line in f:
            t = line.rstrip("\n").split("\t")
            try: snp = t[idx["snp"]]
            except IndexError: continue
            if snp not in stroke_snps: continue
            try: beta=float(t[idx["beta"]]); se=float(t[idx["se"]])
            except (ValueError,IndexError): continue
            eaf = t[idx["eaf"]] if "eaf" in idx else "NA"
            N = t[idx["n"]] if "n" in idx else (str(n_const) if n_const else "NA")
            rows[snp] = (t[idx["ea"]], t[idx["oa"]], eaf, beta, se, t[idx["p"]], N)
    print(f"[outcome {out}] matched {len(rows)}/{len(stroke_snps)}")
    write_outcome("Stroke", out, rows, out)

# ---- 2. SBP as OUTCOME (chr:pos) for CAD, HF, Stroke exposures ----
SBP_EXPS = ["CAD","HF","Stroke"]
need_rsids = set()
for e in SBP_EXPS: need_rsids.update(clumped_snps(e))
# rsid -> (chr,pos) from EUR bim; and (chr,pos)->rsid
rs2pos = {}
with io.open(BIM, encoding="utf-8", errors="replace") as f:
    for line in f:
        p = line.split()
        if len(p) >= 4 and p[1] in need_rsids:
            rs2pos[p[1]] = (p[0], p[3])
pos2rs = {v: k for k, v in rs2pos.items()}
print(f"[SBP outcome] resolved {len(rs2pos)}/{len(need_rsids)} exposure rsIDs to chr:pos")
sbp_rows = {}
with gzip.open(SBP_RAW, "rt", errors="replace") as f:
    hdr = f.readline().split()
    ix = {c: hdr.index(c) for c in ["MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","P"]}
    for line in f:
        t = line.split()
        try:
            mk = t[ix["MarkerName"]].split(":"); key = (mk[0], mk[1])
        except IndexError: continue
        if key not in pos2rs: continue
        rs = pos2rs[key]
        try: beta=float(t[ix["Effect"]]); se=float(t[ix["StdErr"]])
        except (ValueError,IndexError): continue
        sbp_rows[rs] = (t[ix["Allele1"]].upper(), t[ix["Allele2"]].upper(), t[ix["Freq1"]],
                        beta, se, t[ix["P"]], str(SBP_N))
print(f"[SBP outcome] matched {len(sbp_rows)} SNPs in SBP GWAS")
for e in SBP_EXPS:
    write_outcome(e, "SBP", sbp_rows, "SBP")
print("done")
