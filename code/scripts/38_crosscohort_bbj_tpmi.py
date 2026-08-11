# -*- coding: utf-8 -*-
"""Step 38: CROSS-COHORT East Asian MR — Biobank Japan exposure -> TPMI (Taiwan) outcome.

Peer review (C5) noted the EAS analyses are within-cohort (BBJ->BBJ, TPMI->TPMI), i.e. overlapping-
sample MR that can bias toward confounded observational effects. This builds the zero-overlap
alternative: instruments from BBJ (Japanese) tested against outcomes from TPMI (Taiwanese) — two
populations with no shared participants, so it is a clean two-sample design for the headline edges.

Join: BBJ instruments are hg19 rsID; TPMI is hg38 chr:pos (MarkerID=chr_pos_A1_A2, effect=Allele2,
no rsID). Bridge = rsID -> hg19 pos (1000G-EAS bim) -> liftover hg19->hg38 -> match TPMI MarkerID at
that position with allele check. Outcome files are written keyed by the BBJ rsID so 06-style MR
harmonises by rsID. Resumable: per-outcome files skipped if present; the rsID->hg38 map is cached.
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

import gzip, io, os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyliftover import LiftOver

ROOT = CKM_ROOT
BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
UP = os.path.join(BASE, "data/uploads/tpmi")
BIM = os.path.join(BASE, "data/g1000_eas_merged.bim")     # hg19
HARM = os.path.join(BASE, "data/harmonised_crosscohort"); os.makedirs(HARM, exist_ok=True)
MAPCACHE = os.path.join(BASE, "data/bbj_rsid_hg38.pkl")

def log(m): print(m, flush=True)

CROSS_EXPS = ["BMI", "SBP", "LDL", "HDL", "TC", "TG", "HbA1c", "T2D", "CAD"]   # BBJ .eas.clumped
CROSS_OUTS = {"CAD": "411.4", "HF": "428", "Stroke": "433.21", "T2D": "250.2", "CKD": "585.3"}  # TPMI PheCodes

def bbj_instruments(node):
    f = os.path.join(INSTR, f"{node}.eas.clumped.tsv")
    if not os.path.exists(f):
        return {}
    out = {}
    with io.open(f, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: hdr.index(c) for c in ("SNP", "effect_allele", "other_allele")}
        for ln in fh:
            t = ln.rstrip("\n").split("\t")
            out[t[ix["SNP"]]] = (t[ix["effect_allele"]].upper(), t[ix["other_allele"]].upper())
    return out

# ---- 1. rsID -> hg38 (chr,pos) for the union of BBJ instruments (cached) ----
def build_rsid_hg38():
    if os.path.exists(MAPCACHE):
        return pickle.load(open(MAPCACHE, "rb"))
    union = set()
    for e in CROSS_EXPS:
        union |= set(bbj_instruments(e))
    log(f"[map] {len(union)} unique BBJ instrument rsIDs")
    # rsID -> hg19 (chr,pos) via EAS bim
    hg19 = {}
    with io.open(BIM, encoding="utf-8", errors="replace") as f:
        for ln in f:
            p = ln.split("\t") if "\t" in ln else ln.split()
            if len(p) < 4:
                continue
            if p[1] in union and p[1] not in hg19:
                hg19[p[1]] = (p[0], int(p[3]))
    log(f"[map] {len(hg19)}/{len(union)} mapped to hg19 via EAS bim")
    # liftover hg19 -> hg38
    lo = LiftOver("hg19", "hg38")
    hg38 = {}
    for rs, (c, pos) in hg19.items():
        r = lo.convert_coordinate("chr" + c, pos)
        if r:
            hg38[rs] = (c.replace("chr", ""), r[0][1])   # (chr, hg38pos)
    log(f"[map] {len(hg38)}/{len(hg19)} lifted hg19->hg38")
    pickle.dump(hg38, open(MAPCACHE, "wb"))
    return hg38

HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"

def comp(a):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(a, a)

def alleles_match(bea, boa, a1, a2):
    s1, s2 = {bea, boa}, {a1, a2}
    return s1 == s2 or s1 == {comp(a1), comp(a2)}

def main():
    hg38 = build_rsid_hg38()
    # (chr,pos) -> list of rsIDs wanting it
    want = {}
    for rs, cp in hg38.items():
        want.setdefault(cp, []).append(rs)

    for out, phecode in CROSS_OUTS.items():
        # one outcome file PER exposure, but we stream the TPMI outcome ONCE for all rsIDs
        exp_files = {e: os.path.join(HARM, f"BBJ_{e}__TPMI_{out}.outcome.tsv") for e in CROSS_EXPS}
        if all(os.path.exists(p) for p in exp_files.values()):
            log(f"  *->TPMI_{out}: all present, skip"); continue
        src = os.path.join(UP, f"{phecode}.saige.txt.gz")
        if not os.path.exists(src):
            log(f"  TPMI outcome {out} ({phecode}) not found at {src}"); continue
        # rsID -> (ea, oa, eaf, beta, se, p) from the TPMI outcome at the lifted position
        got = {}
        with gzip.open(src, "rt", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {c: hdr.index(c) for c in
                  ("CHR", "POS", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "p.value", "QC")}
            for line in f:
                t = line.rstrip("\n").split("\t")
                try:
                    cp = (t[ix["CHR"]].replace("chr", ""), int(t[ix["POS"]]))
                except (ValueError, IndexError):
                    continue
                rsids = want.get(cp)
                if not rsids:
                    continue
                # QC predicate. This script DROPS QC=="FAIL"; scripts 20/22/25 KEEP only QC=="PASS".
                # Those differ only if the QC column ever holds a third value. Verified empirically
                # 2026-07-23 by streaming the full QC column of six TPMI SAIGE files
                # (585.3, 250.2, 411.4, 428, 433.21 binary + BMI quantitative; 15,588,353 rows each,
                # 93,530,118 rows total): the ONLY distinct values are "PASS" and "FAIL" — exactly,
                # no case or whitespace variants — and PASS+FAIL sums to the row count in every file.
                # The two predicates are therefore EQUIVALENT on the shipped data; 0 variants differ.
                # If a future TPMI release introduces a third QC value they would diverge, so keep
                # this comment with the code and re-run the tally on any new upload.
                if t[ix["QC"]].strip().upper() == "FAIL":
                    continue
                a1, a2 = t[ix["Allele1"]].upper(), t[ix["Allele2"]].upper()   # effect = Allele2
                try:
                    beta = float(t[ix["BETA"]]); se = float(t[ix["SE"]])
                except (ValueError, IndexError):
                    continue
                for rs in rsids:
                    if rs in got:
                        continue
                    bea, boa = None, None
                    for e in CROSS_EXPS:
                        inst = bbj_instruments(e)
                        if rs in inst:
                            bea, boa = inst[rs]; break
                    if bea and not alleles_match(bea, boa, a1, a2):
                        continue     # position hit but alleles incompatible (different variant)
                    got[rs] = (a2, a1, t[ix["AF_Allele2"]], beta, se, t[ix["p.value"]])
        # write per-exposure files
        for e in CROSS_EXPS:
            p = exp_files[e]
            if os.path.exists(p):
                continue
            inst = bbj_instruments(e)
            m = 0
            with io.open(p, "w", encoding="utf-8") as w:
                w.write(HDR)
                for rs in inst:
                    if rs in got:
                        ea, oa, eaf, b, se, pv = got[rs]
                        w.write(f"{rs}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{pv}\t\tTPMI_{out}\n"); m += 1
            log(f"  BBJ_{e}->TPMI_{out}: {m}/{len(inst)} SNPs")

    n = sum(os.path.exists(os.path.join(HARM, f"BBJ_{e}__TPMI_{o}.outcome.tsv"))
            for e in CROSS_EXPS for o in CROSS_OUTS)
    log(f"\n[done] {n} cross-cohort outcome files (expected {len(CROSS_EXPS) * len(CROSS_OUTS)})")

if __name__ == "__main__":
    main()
