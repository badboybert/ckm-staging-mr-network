# -*- coding: utf-8 -*-
"""Step 30 (ITEM 3 robustness): TWO-SAMPLE population EAS test = KoGES-DM instruments -> TWB-Waist outcome.
No sample overlap (Korea vs Taiwan), so immune to the one-sample bias in the within-KoGES DM->BMI edge.
KoGES = hg19 rsID; TWB = hg38 chr:pos (Regenie, effect allele = ALLELE1, p=10^-LOG10P).
Liftover the DM instruments hg19->hg38, match TWB by chr:pos, write a standard outcome TSV keyed by rsID."""
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
from pyliftover import LiftOver

BASE = P4_BASE
DMRAW = os.path.join(BASE, "data/uploads/koges/phenocode-KoGES_DM.tsv.gz")
TWB = os.path.join(BASE, "data/uploads/twb/twb_irnt_gwas_waist_irnt.regenie")
CLUMP = os.path.join(BASE, "data/instruments/KoGES_DM.eas.clumped.tsv")
DEST = os.path.join(BASE, "data/harmonised_eas/KoGES_DM__TWB_Waist.outcome.tsv")

# 1. instrument rsIDs
want = set()
with io.open(CLUMP, encoding="utf-8") as f:
    next(f)
    for line in f:
        want.add(line.split("\t")[0])

# 2. hg19 chr:pos for each instrument from KoGES DM raw
pos19 = {}   # rsid -> (chrom, pos)
with gzip.open(DMRAW, "rt", errors="replace") as f:
    hdr = f.readline().rstrip("\n").split("\t"); ix = {c: hdr.index(c) for c in hdr}
    for line in f:
        t = line.rstrip("\n").split("\t")
        rs = t[ix["rsids"]].strip().split(",")[0]
        if rs in want:
            pos19[rs] = (t[ix["chrom"]], int(t[ix["pos"]]))
        if len(pos19) == len(want):
            break

# 3. liftover hg19->hg38
lo = LiftOver("hg19", "hg38")
pos38 = {}   # (chrom,pos38) -> rsid
lifted = {}
for rs, (c, p) in pos19.items():
    r = lo.convert_coordinate(f"chr{c}", p - 1)   # pyliftover 0-based
    if r:
        c38 = r[0][0].replace("chr", ""); p38 = r[0][1] + 1
        pos38[(c38, p38)] = rs
        lifted[rs] = (c38, p38)
print(f"instruments: {len(want)} | hg19 found: {len(pos19)} | lifted to hg38: {len(lifted)}")

# 4. scan TWB Waist, match by hg38 chr:pos
twb = {}   # rsid -> dict
with io.open(TWB, encoding="utf-8", errors="replace") as f:
    hdr = f.readline().split(); ix = {c: hdr.index(c) for c in hdr}
    for line in f:
        t = line.split()
        key = (t[ix["CHROM"]], int(t[ix["GENPOS"]]))
        if key in pos38:
            rs = pos38[key]
            try:
                beta = float(t[ix["BETA"]]); se = float(t[ix["SE"]])
                p = 10 ** (-float(t[ix["LOG10P"]]))
            except (ValueError, IndexError):
                continue
            twb[rs] = dict(a1=t[ix["ALLELE1"]], a0=t[ix["ALLELE0"]], af=t[ix["A1FREQ"]],
                           beta=beta, se=se, p=p, n=t[ix["N"]])

print(f"matched in TWB Waist (by hg38 pos): {len(twb)}/{len(lifted)}")
missing = sorted(set(lifted) - set(twb))
if missing:
    print("  not in TWB:", missing)

# 5. write outcome TSV keyed by rsID (effect allele = TWB ALLELE1)
os.makedirs(os.path.dirname(DEST), exist_ok=True)
with io.open(DEST, "w", encoding="utf-8") as w:
    w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
    for rs, d in twb.items():
        w.write(f"{rs}\t{d['a1']}\t{d['a0']}\t{d['af']}\t{d['beta']}\t{d['se']}\t{d['p']}\t{d['n']}\tTWB_Waist\n")
print(f"wrote {DEST}")
