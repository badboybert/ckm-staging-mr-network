# -*- coding: utf-8 -*-
"""Step 43: extract RF -> HF-subtype outcome edges (Henry 2024 HERMES non-ischemic HF subtypes).

Tests the obesity-cardiomyopathy route the paper proposes: does adiposity reach non-ischemic HFpEF
and HFrEF directly, and how does that compare with all-cause HF? Also T2D (diabetic cardiomyopathy)
and CAD (should be attenuated into NON-ischemic subtypes vs all-cause HF, since ischemia is excluded).

Outcomes are stage-4 disease GWAS in lib_formats (niHF, niHFpEF, niHFrEF). Reuses the 05-style
rsID extractor. Run: python 43_hf_subtypes.py  ->  data/harmonised_hfsub/<exp>__<sub>.outcome.tsv
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

import io, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_formats import FORMATS, ROOT, _open

BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_hfsub"); os.makedirs(HARM, exist_ok=True)

EXPS = ["BMI", "SBP", "LDL", "TC", "TG", "HbA1c", "T2D", "CAD"]
SUBS = ["niHF", "niHFpEF", "niHFrEF"]
HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"

def clumped(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f):
        return set()
    with io.open(f, encoding="utf-8") as fh:
        fh.readline()
        return {ln.split("\t", 1)[0] for ln in fh}

need = set()
for e in EXPS:
    need |= clumped(e)

for sub in SUBS:
    if all(os.path.exists(os.path.join(HARM, f"{e}__{sub}.outcome.tsv")) for e in EXPS):
        print(f"  *->{sub}: all present, skip", flush=True); continue
    cfg = FORMATS[sub]; C = cfg["cols"]; path = os.path.join(ROOT, cfg["path"])
    rows = {}
    with _open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: hdr.index(v) for k, v in C.items() if v and v in hdr}
        for line in f:
            t = line.rstrip("\n").split("\t")
            try:
                snp = t[idx["snp"]]
            except IndexError:
                continue
            if snp not in need:
                continue
            try:
                beta = float(t[idx["beta"]]); se = float(t[idx["se"]])
            except (ValueError, IndexError):
                continue
            eaf = t[idx["eaf"]] if "eaf" in idx else "NA"
            N = t[idx["n"]] if "n" in idx else "NA"
            rows[snp] = (t[idx["ea"]], t[idx["oa"]], eaf, beta, se, t[idx["p"]], N)
    print(f"[{sub}] matched {len(rows)}/{len(need)} needed SNPs", flush=True)
    for e in EXPS:
        outf = os.path.join(HARM, f"{e}__{sub}.outcome.tsv")
        if os.path.exists(outf):
            continue
        es = clumped(e); m = 0
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for s in es:
                if s in rows:
                    ea, oa, eaf, b, se, p, N = rows[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{sub}\n"); m += 1
        print(f"  {e}->{sub}: {m}/{len(es)} SNPs", flush=True)
print("done", flush=True)
