# -*- coding: utf-8 -*-
"""Step 70: unadjusted waist-to-hip ratio (Pulit 2019 GIANT+UKBB, EUR, build37) as an adiposity
exposure, for the round-2 MAJOR-7 / Tier-2 #10 triangulation.

Reviewer's point: WHRadjBMI is a BMI-residualised exposure, so genetic adjustment for BMI can induce
collider effects and changes the estimand. A null WHRadjBMI->HF does NOT prove central fat "spares"
HF or that overall mass is the causal component. The recommended triangulation is UNADJUSTED WHR,
which retains the overall-mass component, alongside WHRadjBMI (distribution only) and BMI (mass).

Same release and format as the shipped WHRadjBMI file (CHR POS SNP[rsID:EA:OA] ...), build37,
rsID-keyed after stripping the SNP suffix.

Outcomes tested: CAD, all-cause HF, T2D (rsID-keyed, via lib_formats) AND the all-aetiology
HFpEF/HFrEF from Enzan 2025 (position-keyed, mapped via the build37 bim, as in step 68).

Run:  python 70_whr_unadj.py sig     (Part A -> WHR_unadj.sig.tsv)
      Rscript 04_clump.R WHR_unadj   (clump)
      python 70_whr_unadj.py out     (Part B -> harmonised_whr/WHR_unadj__<out>.outcome.tsv)
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

import gzip, io, os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_formats import FORMATS, ROOT, _open

BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_whr"); os.makedirs(HARM, exist_ok=True)
WHR = os.path.join(BASE, "data/raw/WHR_unadj_Pulit2019_EUR.txt.gz")
BIM = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur.bim")
HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"

RSID_OUTS = ["CAD", "HF", "T2D"]                    # via lib_formats FORMATS
POS_OUTS = {                                        # Enzan 2025, position-keyed build37
    "HFpEF": "data/raw/HF_HFpEF_allcause_GCST90654629_transancestry.tsv",
    "HFrEF": "data/raw/HF_HFrEF_allcause_GCST90654628_transancestry.tsv",
}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def rsid(s):
    return s.split(":", 1)[0]


def part_a_sig(pthr=5e-8):
    out = os.path.join(INSTR, "WHR_unadj.sig.tsv"); k = 0
    with gzip.open(WHR, "rt", errors="replace") as f, io.open(out, "w", encoding="utf-8") as w:
        hdr = f.readline().split()
        ix = {c: hdr.index(c) for c in ("SNP", "Tested_Allele", "Other_Allele",
              "Freq_Tested_Allele", "BETA", "SE", "P", "N")}
        w.write(HDR)
        for line in f:
            t = line.split()
            try:
                p = float(t[ix["P"]])
            except (ValueError, IndexError):
                continue
            if not (p < pthr):
                continue
            rs = rsid(t[ix["SNP"]])
            if not rs.startswith("rs"):
                continue
            try:
                b = float(t[ix["BETA"]]); se = float(t[ix["SE"]])
            except (ValueError, IndexError):
                continue
            w.write(f"{rs}\t{t[ix['Tested_Allele']].upper()}\t{t[ix['Other_Allele']].upper()}\t"
                    f"{t[ix['Freq_Tested_Allele']]}\t{b}\t{se}\t{p}\t{t[ix['N']]}\tWHR_unadj\n"); k += 1
    print(f"[WHR_unadj] wrote {k} genome-wide-sig instruments -> WHR_unadj.sig.tsv", flush=True)


def _clumped():
    cf = os.path.join(INSTR, "WHR_unadj.clumped.tsv")
    if not os.path.exists(cf):
        sys.exit("run: Rscript 04_clump.R WHR_unadj  first")
    out = {}
    with io.open(cf, encoding="utf-8") as fh:
        fh.readline()
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            out[p[0]] = (p[1].upper(), p[2].upper())   # rsID -> (EA, OA)
    return out


def _alleles_match(e_ea, e_oa, o_ea, o_oa):
    if not all(a in COMP for a in (e_ea, e_oa, o_ea, o_oa)):
        return False
    e, o = {e_ea, e_oa}, {o_ea, o_oa}
    return e == o or {COMP[e_ea], COMP[e_oa]} == o


def part_b_outcomes():
    snps = _clumped()
    print(f"[WHR_unadj] {len(snps)} clumped instruments", flush=True)

    # rsID-keyed outcomes via lib_formats
    for o in RSID_OUTS:
        outf = os.path.join(HARM, f"WHR_unadj__{o}.outcome.tsv")
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
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for s in snps:
                if s in rows:
                    ea, oa, eaf, b, se, p, N = rows[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{o}\n")
        print(f"  WHR_unadj->{o}: {len(rows)}/{len(snps)} SNPs", flush=True)

    # position-keyed all-aetiology subtypes (Enzan): rsID -> build37 pos -> outcome, allele-checked
    pos_of = {}
    with io.open(BIM, encoding="utf-8", errors="replace") as f:
        for ln in f:
            p = ln.split()
            if p[1] in snps:
                pos_of.setdefault(p[1], (p[0], p[3]))
    wanted = {pos_of[rs]: rs for rs in pos_of}
    for o, relpath in POS_OUTS.items():
        rows = {}
        with io.open(os.path.join(BASE, relpath), encoding="utf-8", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {c: hdr.index(c) for c in ["chromosome", "base_pair_location", "effect_allele",
                  "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "n"]}
            for line in f:
                t = line.rstrip("\n").split("\t")
                try:
                    key = (t[ix["chromosome"]], t[ix["base_pair_location"]])
                except IndexError:
                    continue
                rs = wanted.get(key)
                if rs is None or rs in rows:
                    continue
                rows[rs] = (t[ix["effect_allele"]].upper(), t[ix["other_allele"]].upper(),
                            t[ix["effect_allele_frequency"]], t[ix["beta"]],
                            t[ix["standard_error"]], t[ix["p_value"]], t[ix["n"]])
        outf = os.path.join(HARM, f"WHR_unadj__{o}.outcome.tsv")
        kept = 0
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write(HDR)
            for s, (e_ea, e_oa) in snps.items():
                if s not in rows:
                    continue
                o_ea, o_oa, eaf, b, se, p, N = rows[s]
                if not _alleles_match(e_ea, e_oa, o_ea, o_oa):
                    continue
                w.write(f"{s}\t{o_ea}\t{o_oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{o}\n"); kept += 1
        print(f"  WHR_unadj->{o}: {kept}/{len(snps)} SNPs (position-mapped, allele-checked)", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sig"
    (part_a_sig if mode == "sig" else part_b_outcomes)()
