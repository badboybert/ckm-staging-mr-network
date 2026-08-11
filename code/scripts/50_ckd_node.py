# -*- coding: utf-8 -*-
"""Step 50 (2nd tranche): add binary CKD (Wuttke 2019 CKDGen, GCST008065, EUR, build37) as a
kidney-damage node distinct from the quantitative eGFR node (reviewer C2).

CKD is SPACE-delimited: Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum
(effect allele = Allele1, log-OR). rsID-keyed, so no position mapping.

Roles: X->CKD (CKD as OUTCOME, RF -> kidney damage) and CKD->Y (CKD as EXPOSURE -> CVD). CKD is a
STAGE-2 node (kidney disease is a stage-2 CKM component alongside metabolic risk factors and eGFR).
Part A writes CKD.sig.tsv for clumping; Part B writes the outcome files. Run:
  python 50_ckd_node.py            (both: X->CKD outcomes + CKD.sig)
  Rscript 04_clump.R CKD           (clump CKD instruments)
  python 50_ckd_node.py ckdexp     (CKD->Y outcome files, after clumping)
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

BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_ckd"); os.makedirs(HARM, exist_ok=True)
CKD = os.path.join(ROOT, "paper 6/analysis/data/outcomes/CKD_GCST008065.txt.gz")
HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"

X_EXPS = ["BMI", "SBP", "LDL", "HDL", "TC", "TG", "HbA1c", "T2D", "CAD", "HF"]   # X -> CKD
CKD_OUTS = ["CAD", "HF", "Stroke"]                                              # CKD -> Y

def log(m): print(m, flush=True)

def clumped(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f):
        return set()
    with io.open(f, encoding="utf-8") as fh:
        fh.readline()
        return {ln.split("\t", 1)[0] for ln in fh}

def stream_ckd():
    """yield (rsid, ea, oa, eaf, beta, se, p, n) from the space-delimited CKD GWAS."""
    with gzip.open(CKD, "rt", errors="replace") as f:
        f.readline()
        for line in f:
            t = line.split()
            if len(t) < 10:
                continue
            yield (t[2], t[3].upper(), t[4].upper(), t[5], t[6], t[7], t[8], t[9])

def x_to_ckd():
    need = set()
    for e in X_EXPS:
        need |= clumped(e)
    rows = {}
    for rs, ea, oa, eaf, b, se, p, n in stream_ckd():
        if rs in need:
            rows[rs] = (ea, oa, eaf, b, se, p, n)
    log(f"[CKD outcome] matched {len(rows)}/{len(need)} instrument rsIDs")
    for e in X_EXPS:
        outf = os.path.join(HARM, f"{e}__CKD.outcome.tsv")
        es = clumped(e); m = 0
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for s in es:
                if s in rows:
                    ea, oa, eaf, b, se, p, n = rows[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{n}\tCKD\n"); m += 1
        log(f"  {e}->CKD: {m}/{len(es)} SNPs")

def ckd_sig(pthr=5e-8):
    out = os.path.join(INSTR, "CKD.sig.tsv"); k = 0
    with io.open(out, "w", encoding="utf-8") as w:
        w.write(HDR)
        for rs, ea, oa, eaf, b, se, p, n in stream_ckd():
            try:
                if float(p) >= pthr:
                    continue
            except ValueError:
                continue
            if not rs.startswith("rs"):
                continue
            w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{n}\tCKD\n"); k += 1
    log(f"[CKD exposure] wrote {k} genome-wide-sig instruments -> CKD.sig.tsv (clump: Rscript 04_clump.R CKD)")

def ckd_to_y():
    cf = os.path.join(INSTR, "CKD.clumped.tsv")
    if not os.path.exists(cf):
        sys.exit("clump CKD first: Rscript 04_clump.R CKD")
    with io.open(cf, encoding="utf-8") as fh:
        fh.readline()
        snps = {ln.split("\t", 1)[0] for ln in fh}
    log(f"[CKD->Y] {len(snps)} clumped CKD instruments")
    for o in CKD_OUTS:
        cfg = FORMATS[o]; C = cfg["cols"]; path = os.path.join(ROOT, cfg["path"]); nconst = cfg.get("n_const")
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
                    b = float(t[idx["beta"]]); se = float(t[idx["se"]])
                except (ValueError, IndexError):
                    continue
                eaf = t[idx["eaf"]] if "eaf" in idx else "NA"
                N = t[idx["n"]] if "n" in idx else (str(nconst) if nconst else "NA")
                rows[snp] = (t[idx["ea"]], t[idx["oa"]], eaf, b, se, t[idx["p"]], N)
        outf = os.path.join(HARM, f"CKD__{o}.outcome.tsv"); m = 0
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for s in snps:
                if s in rows:
                    ea, oa, eaf, b, se, p, N = rows[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{o}\n"); m += 1
        log(f"  CKD->{o}: {m}/{len(snps)} SNPs")

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "out"
    if mode == "ckdexp":
        ckd_to_y()
    else:
        x_to_ckd(); ckd_sig()
