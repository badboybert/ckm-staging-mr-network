# -*- coding: utf-8 -*-
"""SBP (Evangelou) has no rsID (MarkerName = chr:pos:type). Extract genome-wide-sig
SNPs, then map (chr,pos)->rsid via the hg19 1000G EUR .bim (both build37)."""
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
SBP  = os.path.join(ROOT, "paper 4/independent_build/data/raw/SBP_Evangelou2018_EUR.txt.gz")
BIM  = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur.bim")
OUT  = os.path.join(ROOT, "paper 4/independent_build/data/instruments/SBP.sig.tsv")
N_SBP = 757601

# Pass 1: collect sig rows keyed by (chr,pos)
sig = {}
with gzip.open(SBP, "rt", errors="replace") as f:
    hdr = f.readline().split()
    ix = {c: hdr.index(c) for c in ["MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","P","TotalSampleSize"]}
    for line in f:
        t = line.split()
        try:
            p = float(t[ix["P"]])
        except (ValueError, IndexError):
            continue
        if not (p < 5e-8): continue
        mk = t[ix["MarkerName"]].split(":")
        if len(mk) < 2: continue
        chrom, pos = mk[0], mk[1]
        try:
            beta = float(t[ix["Effect"]]); se = float(t[ix["StdErr"]])
        except (ValueError, IndexError):
            continue
        sig[(chrom, pos)] = (t[ix["Allele1"]].upper(), t[ix["Allele2"]].upper(),
                             t[ix["Freq1"]], beta, se, p)
print(f"[SBP] {len(sig)} genome-wide-sig positions", flush=True)

# Pass 2: map (chr,pos)->rsid via bim
rsid_map = {}
with io.open(BIM, encoding="utf-8", errors="replace") as f:
    for line in f:
        c, rs, cm, bp, a1, a2 = line.split("\t") if "\t" in line else line.split()
        key = (c, bp)
        if key in sig and key not in rsid_map:
            rsid_map[key] = rs
print(f"[SBP] mapped {len(rsid_map)}/{len(sig)} to rsIDs via bim", flush=True)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with io.open(OUT, "w", encoding="utf-8") as w:
    w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
    k = 0
    for key, rs in rsid_map.items():
        ea, oa, eaf, beta, se, p = sig[key]
        w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{beta}\t{se}\t{p}\t{N_SBP}\tSBP\n"); k += 1
print(f"[SBP] wrote {k} instruments -> {os.path.basename(OUT)}", flush=True)
