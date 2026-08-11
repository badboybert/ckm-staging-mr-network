# -*- coding: utf-8 -*-
"""Step 40: add coronary-artery-calcium (CAC; Kavousi 2023, GCST90278456, EUR, GRCh37) as a
STAGE-3 (subclinical CVD) node to the European network.

Peer review C2 / the 2026 CKM guideline: the network jumps stage 2 (metabolic) -> stage 4 (clinical
CVD) with no stage-3 subclinical tier. CAC is the guideline's named stage-3 marker. This lets the
staging concordance test span 1 -> 2 -> 3 -> 4.

CAC is GRCh37 chr:pos with no rsID (variant_id = chr_pos_ea_oa), same build as the EUR network, so no
liftover. Two roles:
  X -> CAC  (CAC as OUTCOME): map each exposure instrument rsID -> (chr,pos) via the build37 EUR bim,
            read CAC at that position. Robust (many instruments).
  CAC -> Y  (CAC as EXPOSURE): extract CAC genome-wide-sig variants, map (chr,pos) -> rsID, clump
            (04_clump.R), then read each disease outcome. Thin (~few loci) but gives the stage 3->4 link.

Resumable: outcome files skipped if present; the sig-extraction writes CAC.sig.tsv for clumping.
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

import gzip, io, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = CKM_ROOT
BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised")
BIM = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur.bim")           # build37
CAC = os.path.join(BASE, "data/raw/CAC_Kavousi2023_EUR.tsv.gz")

def log(m): print(m, flush=True)

# exposures whose instruments we read from CAC (CAC as outcome). Stage 1/2 forward + disease for completeness.
OUT_EXPS = ["BMI", "SBP", "LDL", "HDL", "TC", "TG", "HbA1c", "T2D", "eGFR", "CAD", "HF", "Stroke"]
HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"

def clumped_rsids(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f):
        return set()
    with io.open(f, encoding="utf-8") as fh:
        fh.readline()
        return {ln.split("\t", 1)[0] for ln in fh}

def read_cac():
    """CAC keyed by (chr,pos) -> (ea,oa,eaf,beta,se,p,n); plus the set of genome-wide-sig (chr,pos)."""
    at, sig = {}, set()
    with gzip.open(CAC, "rt", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {c: hdr.index(c) for c in ("chromosome", "base_pair_location", "effect_allele",
              "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "n")}
        for line in f:
            t = line.rstrip("\n").split("\t")
            try:
                key = (t[ix["chromosome"]], t[ix["base_pair_location"]])
                p = float(t[ix["p_value"]])
            except (ValueError, IndexError):
                continue
            at[key] = (t[ix["effect_allele"]].upper(), t[ix["other_allele"]].upper(),
                       t[ix["effect_allele_frequency"]], t[ix["beta"]], t[ix["standard_error"]],
                       t[ix["p_value"]], t[ix["n"]])
            if p < 5e-8:
                sig.add(key)
    log(f"[CAC] {len(at)} variants, {len(sig)} genome-wide-sig positions")
    return at, sig

def main():
    at, sig = read_cac()
    # exposure instruments we need to place (rsID) + CAC-sig positions we need to name (chr,pos)
    exp_rsids = set()
    for e in OUT_EXPS:
        exp_rsids |= clumped_rsids(e)
    # one bim pass: rsID->(chr,pos) for exposure instruments; (chr,pos)->rsID for CAC-sig
    rs_to_cp, cp_to_rs = {}, {}
    with io.open(BIM, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.split("\t") if "\t" in line else line.split()
            if len(p) < 4:
                continue
            rs, chrom, pos = p[1], p[0], p[3]
            if rs in exp_rsids and rs not in rs_to_cp:
                rs_to_cp[rs] = (chrom, pos)
            if (chrom, pos) in sig and (chrom, pos) not in cp_to_rs:
                cp_to_rs[(chrom, pos)] = rs
    log(f"[bim] {len(rs_to_cp)}/{len(exp_rsids)} exposure instruments placed; "
        f"{len(cp_to_rs)}/{len(sig)} CAC-sig positions named")

    # ---- Part A: CAC as EXPOSURE — write CAC.sig.tsv (for 04_clump.R) ----
    sigf = os.path.join(INSTR, "CAC.sig.tsv")
    k = 0
    with io.open(sigf, "w", encoding="utf-8") as w:
        w.write(HDR)
        for cp, rs in cp_to_rs.items():
            ea, oa, eaf, b, se, p, n = at[cp]
            w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{n}\tCAC\n"); k += 1
    log(f"[CAC exposure] wrote {k} genome-wide-sig instruments -> CAC.sig.tsv (clump next: Rscript 04_clump.R CAC)")

    # ---- Part B: CAC as OUTCOME — X -> CAC ----
    for e in OUT_EXPS:
        outf = os.path.join(HARM, f"{e}__CAC.outcome.tsv")
        if os.path.exists(outf):
            log(f"  {e}->CAC: exists, skip"); continue
        rsids = clumped_rsids(e)
        m = 0
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for rs in rsids:
                cp = rs_to_cp.get(rs)
                if cp and cp in at:
                    ea, oa, eaf, b, se, p, n = at[cp]
                    w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{n}\tCAC\n"); m += 1
        log(f"  {e}->CAC: {m}/{len(rsids)} SNPs -> {os.path.basename(outf)}")

if __name__ == "__main__":
    main()
