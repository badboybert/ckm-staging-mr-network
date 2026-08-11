# -*- coding: utf-8 -*-
"""Step 37: generate the 21 directed edges the original network never estimated, so the graph
covers all 132 non-self directions among the 12 traits instead of 111.

The gap has one mechanical cause: 05_extract_outcomes.py listed SBP as an exposure only and Stroke
as an outcome only. So no X->SBP edge and no Stroke->X edge was ever extracted. This script fills
exactly those, writing outcome files in the identical schema, so 06_mr.R then computes them with the
same harmonisation and estimators as the existing 111.

Two mechanically different cases:
  Stroke -> X : Stroke has 24 clumped instruments already; extract them from each outcome GWAS
                (all in lib_formats.FORMATS).
  X -> SBP    : SBP must be an OUTCOME, but the Evangelou GWAS is keyed by chr:pos, not rsID. Map
                each exposure instrument rsID -> (chr,pos) via the build37 1000G bim (same build as
                Evangelou), then read SBP at those positions.

Resumable: each per-edge outcome file is written once and skipped on re-run; a kill costs at most
one outcome stream. Progress is printed per outcome, unbuffered.
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

import os, io, sys, gzip
sys.path.insert(0, os.path.dirname(__file__))
from lib_formats import FORMATS, ROOT, _open

BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised")
BIM = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur.bim")
SBP_GWAS = os.path.join(BASE, "data/raw/SBP_Evangelou2018_EUR.txt.gz")
SBP_N = 757601
os.makedirs(HARM, exist_ok=True)

def log(m):
    print(m, flush=True)

def load_clumped(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f):
        return {}
    out = {}
    with io.open(f, encoding="utf-8") as fh:
        fh.readline()
        for ln in fh:
            out[ln.split("\t", 1)[0]] = ln  # rsid -> full line (unused, kept for count)
    return out

# ---- the 21 new edges ----
SBP_EXPS = ["BMI", "TG", "HDL", "TC", "LDL", "HbA1c", "CAD", "HF", "T2D", "eGFR", "Stroke"]  # X->SBP
STROKE_OUTS = ["CAD", "HF", "T2D", "eGFR", "BMI", "HDL", "TC", "TG", "LDL", "HbA1c"]          # Stroke->X

HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"

# ============================ Stroke -> X ============================
def do_stroke_edges():
    stroke_snps = set(load_clumped("Stroke"))
    log(f"[Stroke] {len(stroke_snps)} clumped instruments")
    for o in STROKE_OUTS:
        outf = os.path.join(HARM, f"Stroke__{o}.outcome.tsv")
        if os.path.exists(outf):
            log(f"  Stroke->{o}: exists, skip")
            continue
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
                if snp not in stroke_snps:
                    continue
                try:
                    beta = float(t[idx["beta"]]); se = float(t[idx["se"]])
                except (ValueError, IndexError):
                    continue
                eaf = t[idx["eaf"]] if "eaf" in idx else "NA"
                N = t[idx["n"]] if "n" in idx else (str(n_const) if n_const else "NA")
                rows[snp] = (t[idx["ea"]], t[idx["oa"]], eaf, beta, se, t[idx["p"]], N)
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for s in stroke_snps:
                if s in rows:
                    ea, oa, eaf, b, se, p, N = rows[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{o}\n")
        log(f"  Stroke->{o}: {len(rows)}/{len(stroke_snps)} SNPs -> {os.path.basename(outf)}")

# ============================ X -> SBP ============================
def do_sbp_edges():
    # already done?
    if all(os.path.exists(os.path.join(HARM, f"{e}__SBP.outcome.tsv")) for e in SBP_EXPS):
        log("[X->SBP] all present, skip")
        return
    # 1) union of exposure instrument rsIDs
    exp_snps = {e: set(load_clumped(e)) for e in SBP_EXPS}
    union = set().union(*exp_snps.values())
    log(f"[X->SBP] {len(union)} unique exposure instruments across {len(SBP_EXPS)} exposures")

    # 2) rsID -> (chr,pos) from the build37 bim (same build as Evangelou)
    pos_of = {}
    with io.open(BIM, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.split("\t") if "\t" in line else line.split()
            if len(p) < 4:
                continue
            rs = p[1]
            if rs in union and rs not in pos_of:
                pos_of[rs] = (p[0], p[3])
    log(f"[X->SBP] mapped {len(pos_of)}/{len(union)} rsIDs to build37 positions")
    wanted = {}                       # (chr,pos) -> rsid
    for rs, cp in pos_of.items():
        wanted.setdefault(cp, rs)

    # 3) stream SBP once, collect rows at wanted positions
    sbp_at = {}                       # rsid -> (ea,oa,eaf,beta,se,p)
    with gzip.open(SBP_GWAS, "rt", errors="replace") as f:
        hdr = f.readline().split()
        ix = {c: hdr.index(c) for c in
              ["MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P"]}
        for line in f:
            t = line.split()
            try:
                mk = t[ix["MarkerName"]].split(":")
            except IndexError:
                continue
            if len(mk) < 2:
                continue
            key = (mk[0], mk[1])
            rs = wanted.get(key)
            if rs is None or rs in sbp_at:
                continue
            try:
                beta = float(t[ix["Effect"]]); se = float(t[ix["StdErr"]])
            except (ValueError, IndexError):
                continue
            sbp_at[rs] = (t[ix["Allele1"]].upper(), t[ix["Allele2"]].upper(),
                          t[ix["Freq1"]], beta, se, t[ix["P"]])
    log(f"[X->SBP] SBP has data for {len(sbp_at)}/{len(pos_of)} mapped instruments")

    # 4) write per-exposure SBP outcome files
    for e in SBP_EXPS:
        outf = os.path.join(HARM, f"{e}__SBP.outcome.tsv")
        if os.path.exists(outf):
            log(f"  {e}->SBP: exists, skip"); continue
        m = 0
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for s in exp_snps[e]:
                if s in sbp_at:
                    ea, oa, eaf, b, se, p = sbp_at[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{SBP_N}\tSBP\n"); m += 1
        log(f"  {e}->SBP: {m}/{len(exp_snps[e])} SNPs -> {os.path.basename(outf)}")

if __name__ == "__main__":
    do_stroke_edges()
    do_sbp_edges()
    # summary
    n = sum(os.path.exists(os.path.join(HARM, f"{e}__SBP.outcome.tsv")) for e in SBP_EXPS) + \
        sum(os.path.exists(os.path.join(HARM, f"Stroke__{o}.outcome.tsv")) for o in STROKE_OUTS)
    log(f"\n[done] {n} new outcome files present (expected 21)")
