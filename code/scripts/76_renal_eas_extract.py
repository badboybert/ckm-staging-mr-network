# -*- coding: utf-8 -*-
"""Step 76a: harmonised East-Asian CKD outcome files for the two renal contrasts.

WHY THIS EXISTS (round-7 item C018/C102, PI decision Q3 2026-09-16): the SBP->CKD and BMI->CKD
contrasts quoted in Results/Discussion/Supplementary Figure S6a came from 09_zheng_ckd.py, a
self-contained FIXED-effect Python IVW, while every other edge in the paper is TwoSampleMR with
multiplicative random effects (06_mr.R). Two estimators for one claim is one estimator too many.
This step writes the EAS outcome files so 76_renal_eas_mr.R can estimate those two edges with the
SAME battery as the rest of the network; the European side already exists in results/network_ckd.csv
(51_ckd_mr.R) and is NOT recomputed here - there must be exactly one European source.

Outcome: Biobank Japan CKD (Sakaue 2021, GCST90018602; 2,117 cases / 174,345 controls per Table S1).
The file carries no sample-size column, so N is taken from that single declared total and asserted
against Table S1 by the gate, not retyped per row.

Writes: data/harmonised_renal/{SBP,BMI}_eas__CKD.outcome.tsv
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

sys.path.insert(0, os.path.dirname(__file__))

ROOT = CKM_ROOT
BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
OUTDIR = os.path.join(BASE, "data/harmonised_renal")
EAS_CKD = os.path.join(ROOT, "paper 6/analysis/data/outcomes/eas/CKD_GCST90018602.gz")
N_EAS_CKD = 2117 + 174345          # Table S1: Biobank Japan CKD, cases + controls
os.makedirs(OUTDIR, exist_ok=True)


def instrument_rsids(path):
    rs = set()
    with io.open(path, encoding="utf-8") as fh:
        fh.readline()
        for ln in fh:
            t = ln.rstrip("\n").split("\t")
            if t and t[0]:
                rs.add(t[0])
    return rs


def main():
    exps = ["SBP", "BMI"]
    want = {}
    for e in exps:
        want[e] = instrument_rsids(os.path.join(INSTR, f"{e}.eas.clumped.tsv"))
        print(f"{e}: {len(want[e])} EAS instruments", flush=True)
    allrs = set().union(*want.values())

    # chr:pos -> rsid, because the atlas file leaves variant_id NA for most rows
    from lib_formats_eas import BIM
    pos2rs = {}
    with io.open(BIM, errors="replace") as f:
        for ln in f:
            c, rs, cm, bp, a1, a2 = ln.split()
            pos2rs[(c, bp)] = rs
    print(f"bim positions: {len(pos2rs):,}", flush=True)

    found = {}
    with gzip.open(EAS_CKD, "rt", errors="replace") as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = {c: h.index(c) for c in ["chromosome", "base_pair_location", "effect_allele",
                                      "other_allele", "effect_allele_frequency", "beta",
                                      "standard_error", "p_value", "variant_id"]}
        for ln in f:
            t = ln.rstrip("\n").split("\t")
            rs = t[ix["variant_id"]].strip()
            if not rs.startswith("rs"):
                rs = pos2rs.get((t[ix["chromosome"]].strip(), t[ix["base_pair_location"]].strip()))
            if not rs or rs not in allrs or rs in found:
                continue
            try:
                found[rs] = (t[ix["effect_allele"]].upper(), t[ix["other_allele"]].upper(),
                             t[ix["effect_allele_frequency"]], float(t[ix["beta"]]),
                             float(t[ix["standard_error"]]), t[ix["p_value"]])
            except ValueError:
                continue
    print(f"matched {len(found)} of {len(allrs)} instrument variants in the BBJ CKD file", flush=True)

    for e in exps:
        out = os.path.join(OUTDIR, f"{e}_eas__CKD.outcome.tsv")
        n = 0
        with io.open(out, "w", encoding="utf-8", newline="") as fh:
            fh.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
            for rs in sorted(want[e]):
                if rs not in found:
                    continue
                ea, oa, eaf, b, se, p = found[rs]
                fh.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N_EAS_CKD}\tCKD_EAS\n")
                n += 1
        print(f"wrote {out} ({n} variants)", flush=True)
        if n < 3:
            raise SystemExit(f"FAIL: {e} -> CKD has {n} matched variants; MR needs >= 3")


if __name__ == "__main__":
    main()
