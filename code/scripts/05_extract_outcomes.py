# -*- coding: utf-8 -*-
"""Step 5a: for each edge, extract the OUTCOME GWAS rows matching the EXPOSURE's
clumped instrument rsIDs. Streams each outcome file once (grouped by outcome)."""
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

import os, io, gzip, sys
sys.path.insert(0, os.path.dirname(__file__))
from lib_formats import FORMATS, ROOT, _open

BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
HARM  = os.path.join(BASE, "data/harmonised")
os.makedirs(HARM, exist_ok=True)

# FULL BIDIRECTIONAL matrix: exposures = any with .clumped.tsv; outcomes = any clean rsID+beta node
# (adds risk factors as OUTCOMES -> reverse disease->RF edges = the falsification test).
EXPS = ["BMI","SBP","TG","HDL","TC","LDL","HbA1c","CAD","HF","T2D","eGFR"]
OUTS = ["CAD","HF","Stroke","T2D","eGFR","BMI","HDL","TC","TG","LDL","HbA1c"]
EDGES = [(e, o) for e in EXPS for o in OUTS if e != o]

def load_clumped_snps(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f): return set()
    with io.open(f, encoding="utf-8") as fh:
        fh.readline()
        return {ln.split("\t")[0] for ln in fh}

# group needed SNPs by outcome
need = {}
for e, o in EDGES:
    need.setdefault(o, set()).update(load_clumped_snps(e))

# stream each outcome once, collect standardized rows for needed SNPs
for o, snps in need.items():
    cfg = FORMATS[o]; C = cfg["cols"]; path = os.path.join(ROOT, cfg["path"])
    n_const = cfg.get("n_const")
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
            if snp not in snps:
                continue
            try:
                beta = float(t[idx["beta"]]); se = float(t[idx["se"]])
            except (ValueError, IndexError):
                continue
            eaf = t[idx["eaf"]] if "eaf" in idx else "NA"
            N = t[idx["n"]] if "n" in idx else (str(n_const) if n_const else "NA")
            rows[snp] = (t[idx["ea"]], t[idx["oa"]], eaf, beta, se, t[idx["p"]], N)
    print(f"[outcome] {o}: matched {len(rows)}/{len(snps)} needed SNPs", flush=True)
    # write per-edge outcome files
    for e, oo in EDGES:
        if oo != o: continue
        exp_snps = load_clumped_snps(e)
        outf = os.path.join(HARM, f"{e}__{o}.outcome.tsv")
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
            m = 0
            for s in exp_snps:
                if s in rows:
                    ea, oa, eaf, b, se, p, N = rows[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{o}\n"); m += 1
        print(f"   {e}->{o}: {m} outcome SNPs -> {os.path.basename(outf)}", flush=True)
print("done")
