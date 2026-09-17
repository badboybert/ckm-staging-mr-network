# -*- coding: utf-8 -*-
"""Step 79: give the European HbA1c dataset its true sample size.

PROVENANCE CORRECTION (2026-09-16, PI decision: keep the analysed file, correct its description).
scripts/02_download.sh line 16 fetches GWAS Catalog accession GCST90014006 and saves it as
"HbA1c_MAGIC_EUR.h.tsv.gz". GCST90014006 is Mbatchou et al. 2021 (PMID 34017140), the REGENIE
paper's UK Biobank HbA1c GWAS in 389,889 European participants - NOT the MAGIC release of Chen et
al. 2021 (146,806) that the filename asserts. The harmonised file carries no sample-size column, so
lib_formats.py supplied the wrong constant, and every derived file inherited it.

Effect sizes and standard errors are untouched by this: they are whatever the file holds. The one
quantity that depends on N is the Steiger directionality test, which gates the staging graph. This
step rewrites the sample size in the extracted instrument and harmonised outcome files so the gate
runs on the true N; 06_mr.R is then re-run and its output diffed against the previous table.

Writes (in place, after asserting the old value): data/instruments/HbA1c.{clumped,sig}.tsv and
data/harmonised/*__HbA1c.outcome.tsv
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

import glob, io, os, sys

sys.stdout.reconfigure(encoding="utf-8")
BASE = P4_BASE
N_WRONG, N_TRUE = "146806", "389889"


def fix(path):
    with io.open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")
    head = lines[0].split("\t")
    if "N" not in head:
        return 0, "no N column"
    j = head.index("N")
    n_changed = n_other = 0
    for i in range(1, len(lines)):
        if not lines[i].strip():
            continue
        t = lines[i].split("\t")
        if len(t) <= j:
            continue
        if t[j] == N_WRONG:
            t[j] = N_TRUE
            lines[i] = "\t".join(t)
            n_changed += 1
        else:
            n_other += 1
    if n_changed:
        io.open(path, "w", encoding="utf-8", newline="\n").write("\n".join(lines))
    return n_changed, f"{n_other} row(s) carried a different N (left alone)"


def main():
    targets = [os.path.join(BASE, "data/instruments", f"HbA1c.{k}.tsv") for k in ("clumped", "sig")]
    targets += sorted(glob.glob(os.path.join(BASE, "data/harmonised", "*__HbA1c.outcome.tsv")))
    total = 0
    for p in targets:
        if not os.path.exists(p):
            print(f"  MISSING {os.path.basename(p)}")
            continue
        n, note = fix(p)
        total += n
        print(f"  {os.path.basename(p):32s} {n:6d} rows -> N={N_TRUE}   ({note})")
    print(f"\nrewrote {total} sample-size cells across {len(targets)} files")
    # the EAS and TPMI HbA1c files are different datasets and must NOT be touched
    for p in glob.glob(os.path.join(BASE, "data/instruments", "HbA1c.*.tsv")):
        if ".eas." in p or ".tpmi." in p:
            t = io.open(p, encoding="utf-8").read()
            assert N_TRUE not in t, f"{os.path.basename(p)} must not carry the European N"
    print("checked: the East Asian and TPMI HbA1c instrument files are untouched")


if __name__ == "__main__":
    main()
