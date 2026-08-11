# -*- coding: utf-8 -*-
"""Step 28 (ITEM 3): extract genome-wide-sig instruments from KoGES (Korean POPULATION cohort).
KoGES is rsID-keyed, hg19 (verified: rs7903146 @ chr10:114758349), effect allele = alt, beta = log-OR
(binary) or per-SD (continuous, IRNT). Writes standard instrument TSVs so 04_clump_eas.R can clump them.
Exposures for the population-vs-hospital ->BMI triangulation: DM, BMI, WAIST."""
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

BASE = P4_BASE
KDIR = os.path.join(BASE, "data/uploads/koges")
OUT  = os.path.join(BASE, "data/instruments")
os.makedirs(OUT, exist_ok=True)
N_KOGES = 72298   # KoGES sample (continuous N via ac/(2*af) ~ 71.7k; use cohort N)

# node -> (filename, eaf-column-name). continuous uses 'af'; binary DM uses 'control_af'.
KOGES = {
    "KoGES_DM":    ("phenocode-KoGES_DM.tsv.gz",    "control_af"),
    "KoGES_BMI":   ("phenocode-KoGES_BMI.tsv.gz",   "af"),
    "KoGES_WAIST": ("phenocode-KoGES_WAIST.tsv.gz", "af"),
}

def extract(node, pthr=5e-8):
    fn, eaf_col = KOGES[node]
    path = os.path.join(KDIR, fn)
    out = os.path.join(OUT, f"{node}.eas.sig.tsv")
    kept = 0
    with gzip.open(path, "rt", errors="replace") as f, io.open(out, "w", encoding="utf-8") as w:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {c: hdr.index(c) for c in hdr}
        w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
        for line in f:
            t = line.rstrip("\n").split("\t")
            try:
                p = float(t[ix["pval"]])
            except (ValueError, IndexError):
                continue
            if not (p < pthr):
                continue
            try:
                ea = t[ix["alt"]].strip(); oa = t[ix["ref"]].strip()
                beta = float(t[ix["beta"]]); se = float(t[ix["sebeta"]])
            except (ValueError, IndexError):
                continue
            if se <= 0 or math.isnan(beta) or math.isnan(se):
                continue
            rs = t[ix["rsids"]].strip().split(",")[0]
            if not rs.startswith("rs"):
                continue
            eaf = t[ix[eaf_col]] if eaf_col in ix else "NA"
            w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{beta}\t{se}\t{p}\t{N_KOGES}\t{node}\n")
            kept += 1
    print(f"[KoGES extract] {node}: {kept} sig SNPs (p<{pthr:g}) -> {os.path.basename(out)}")
    return kept

if __name__ == "__main__":
    import sys
    nodes = sys.argv[1:] or list(KOGES)
    for n in nodes:
        extract(n)
