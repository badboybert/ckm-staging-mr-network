# -*- coding: utf-8 -*-
"""Step 28b (ITEM 3): build KoGES outcome files for the triangulation edges.
For each edge exp->out, read the clumped exposure rsIDs, scan the full KoGES outcome GWAS,
and write a standard outcome TSV. Within-cohort (same KoGES sample/build) -> harmonise by rsID."""
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
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_eas")
os.makedirs(HARM, exist_ok=True)
N_KOGES = 72298

KFILE = {"KoGES_DM": ("phenocode-KoGES_DM.tsv.gz", "control_af"),
         "KoGES_BMI": ("phenocode-KoGES_BMI.tsv.gz", "af"),
         "KoGES_WAIST": ("phenocode-KoGES_WAIST.tsv.gz", "af")}

EDGES = [("KoGES_DM", "KoGES_BMI"), ("KoGES_DM", "KoGES_WAIST"),
         ("KoGES_BMI", "KoGES_DM"), ("KoGES_WAIST", "KoGES_DM")]

def clumped_rsids(node):
    f = os.path.join(INSTR, f"{node}.eas.clumped.tsv")
    rs = set()
    with io.open(f, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            rs.add(line.split("\t")[0])
    return rs

def build_outcome(exp, out):
    want = clumped_rsids(exp)
    fn, eaf_col = KFILE[out]
    path = os.path.join(KDIR, fn)
    dest = os.path.join(HARM, f"{exp}__{out}.outcome.tsv")
    found = 0
    with gzip.open(path, "rt", errors="replace") as f, io.open(dest, "w", encoding="utf-8") as w:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {c: hdr.index(c) for c in hdr}
        w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
        for line in f:
            t = line.rstrip("\n").split("\t")
            rs = t[ix["rsids"]].strip().split(",")[0]
            if rs not in want:
                continue
            try:
                ea = t[ix["alt"]].strip(); oa = t[ix["ref"]].strip()
                beta = float(t[ix["beta"]]); se = float(t[ix["sebeta"]]); p = float(t[ix["pval"]])
            except (ValueError, IndexError):
                continue
            if se <= 0 or math.isnan(beta):
                continue
            eaf = t[ix[eaf_col]] if eaf_col in ix else "NA"
            w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{beta}\t{se}\t{p}\t{N_KOGES}\t{out}\n")
            found += 1
    print(f"[KoGES outcome] {exp} -> {out}: {found}/{len(want)} instruments found")
    return found

if __name__ == "__main__":
    for exp, out in EDGES:
        build_outcome(exp, out)
