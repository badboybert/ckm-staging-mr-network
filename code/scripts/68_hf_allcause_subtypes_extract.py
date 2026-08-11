# -*- coding: utf-8 -*-
"""Step 68: extract RF -> ALL-AETIOLOGY HFpEF/HFrEF outcome edges (Enzan 2025, PMID 41184235).

Answers round-2 review BLOCKER-2. The shipped subtype analysis uses NON-ISCHEMIC HFpEF/HFrEF (Henry
2025 HERMES, 3,590 / 4,975 cases) as a "specificity control" for a CAD-mediated T2D->HF mechanism -
which the reviewer called structurally circular, because CAD is partly excluded by the non-ischemic
case definition, so "CAD does not reach non-ischemic HF" is built in rather than discovered.

All-aetiology subtypes do not have that problem: CAD is NOT excluded, so a genuinely null CAD->HFpEF
edge (if it occurs) is evidence, not definition. Source = the SAME Enzan 2025 paper the manuscript
already cites for the EAS HF outcome (GCST90668009). all-aetiology, CROSS-ANCESTRY (EUR+EAS) meta-analyses - NOT European-only (only the combined file is
deposited; the "_EUR" label an earlier version used was WRONG):
    GCST90654628 = HFrEF  N=480,269 (23,749 cases: 19,495 EUR + 4,254 EAS)
    GCST90654629 = HFpEF  N=483,263 (26,743 cases: 19,589 EUR + 7,154 EAS)
The EAS component (197,577 EAS controls) OVERLAPS the paper's separately-used Biobank Japan HF outcome
GCST90668009 (same Enzan paper). EUR instruments are applied to a ~42%-EAS outcome, an ancestry/LD
mismatch that ATTENUATES rather than inflates - so the strong CAD/BMI flips below are conservative.
The MTAG variants (GCST90668014/15) are deliberately NOT used - MTAG borrows information across
traits and is not a valid MR outcome.

The Enzan files are keyed by chromosome:base_pair_location with NO rsID (build37, verified: 103/103
BMI instruments position-matched with 103/103 allele concordance). So this maps each instrument
rsID -> (chr,pos) via the build37 1000G bim, then reads the outcome at that position WITH an allele-
concordance check (the SBP-extraction lesson: position-only matching can silently accept a different
variant at the same coordinate).

Run: python 68_hf_allcause_subtypes_extract.py
Out: data/harmonised_hfsub/<exp>__<sub>.outcome.tsv    (sub in {HFpEF, HFrEF})
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

import io, os, sys, collections

BASE = P4_BASE
ROOT = CKM_ROOT
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_hfsub"); os.makedirs(HARM, exist_ok=True)
BIM = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur.bim")

EXPS = ["BMI", "SBP", "LDL", "TC", "TG", "HbA1c", "T2D", "CAD"]
SUBS = {
    "HFpEF": ("data/raw/HF_HFpEF_allcause_GCST90654629_transancestry.tsv", 19589),
    "HFrEF": ("data/raw/HF_HFrEF_allcause_GCST90654628_transancestry.tsv", 19495),
}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"


def clumped(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f):
        return {}
    out = {}
    with io.open(f, encoding="utf-8") as fh:
        fh.readline()
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            out[p[0]] = (p[1].upper(), p[2].upper())   # rsID -> (EA, OA) from the instrument side
    return out


def alleles_match(e_ea, e_oa, o_ea, o_oa):
    """True if the outcome record is the SAME biallelic SNV as the instrument (strand flip allowed)."""
    if not all(a in COMP for a in (e_ea, e_oa, o_ea, o_oa)):
        return False
    e, o = {e_ea, e_oa}, {o_ea, o_oa}
    return e == o or {COMP[e_ea], COMP[e_oa]} == o


def main():
    inst = {e: clumped(e) for e in EXPS}
    union = set().union(*[set(d) for d in inst.values()])
    print(f"union of instruments across {len(EXPS)} exposures: {len(union)}", flush=True)

    # rsID -> (chr,pos); detect position collisions (a coordinate carrying >1 instrument rsID)
    pos_of, at_pos = {}, collections.defaultdict(list)
    with io.open(BIM, encoding="utf-8", errors="replace") as f:
        for ln in f:
            p = ln.split()
            if p[1] in union:
                pos_of.setdefault(p[1], (p[0], p[3]))
                at_pos[(p[0], p[3])].append(p[1])
    collisions = {k: v for k, v in at_pos.items() if len(v) > 1}
    print(f"mapped {len(pos_of)}/{len(union)} instruments to build37 positions; "
          f"position collisions: {len(collisions)}", flush=True)
    wanted = {pos_of[rs]: rs for rs in pos_of}

    for sub, (relpath, _) in SUBS.items():
        path = os.path.join(BASE, relpath)
        rows, multi = {}, collections.Counter()
        with io.open(path, encoding="utf-8", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {c: hdr.index(c) for c in
                  ["chromosome", "base_pair_location", "effect_allele", "other_allele",
                   "beta", "standard_error", "effect_allele_frequency", "p_value", "n"]}
            for line in f:
                t = line.rstrip("\n").split("\t")
                try:
                    key = (t[ix["chromosome"]], t[ix["base_pair_location"]])
                except IndexError:
                    continue
                rs = wanted.get(key)
                if rs is None:
                    continue
                multi[key] += 1
                if rs in rows:                    # first record at this position wins; count the rest
                    continue
                rows[rs] = (t[ix["effect_allele"]].upper(), t[ix["other_allele"]].upper(),
                            t[ix["effect_allele_frequency"]], t[ix["beta"]],
                            t[ix["standard_error"]], t[ix["p_value"]], t[ix["n"]])
        n_multi = sum(1 for k, c in multi.items() if c > 1)
        print(f"[{sub}] matched {len(rows)}/{len(pos_of)} instrument positions "
              f"({n_multi} positions had >1 outcome record, first kept)", flush=True)

        for e in EXPS:
            outf = os.path.join(HARM, f"{e}__{sub}.outcome.tsv")
            m = kept = dropped = 0
            with io.open(outf, "w", encoding="utf-8") as w:
                w.write(HDR)
                for rs, (e_ea, e_oa) in inst[e].items():
                    if rs not in rows:
                        continue
                    m += 1
                    o_ea, o_oa, eaf, b, se, p, N = rows[rs]
                    if not alleles_match(e_ea, e_oa, o_ea, o_oa):
                        dropped += 1
                        continue
                    w.write(f"{rs}\t{o_ea}\t{o_oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{sub}\n"); kept += 1
            print(f"  {e}->{sub}: {kept}/{len(inst[e])} SNPs "
                  f"({dropped} dropped on allele mismatch of {m} position-matched)", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
