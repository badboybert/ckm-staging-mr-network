# -*- coding: utf-8 -*-
"""Step 45: WHRadjBMI (Pulit 2019 GIANT+UKBB, EUR, build37) as an adiposity-DISTRIBUTION exposure.

Path B: BMI captures overall adiposity; WHRadjBMI isolates central fat distribution independent of
BMI. Tests whether the cascade is driven by fat distribution as well as overall mass. Compare
WHRadjBMI->{CAD,HF,T2D,niHFpEF,niHFrEF} against BMI->.

The SNP column is 'rsID:EA:OA' -> strip to rsID. Build37, rsID-keyed, so no liftover.
Part A writes WHRadjBMI.sig.tsv (clump next: Rscript 04_clump.R WHRadjBMI).
Part B extracts WHRadjBMI instruments from each outcome GWAS (needs WHRadjBMI.clumped.tsv first).

⚠ WHRadjBMI includes UK Biobank -> sample overlap with UKB-containing outcome GWAS (CAD, HF). Report
with the overlap caveat; the direction is more robust than the magnitude.
Run: python 45_whradjbmi.py sig   (Part A)  then  python 45_whradjbmi.py out   (Part B, after clump)
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
from lib_formats import FORMATS, ROOT, _open

BASE = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised_whr"); os.makedirs(HARM, exist_ok=True)
WHR = os.path.join(BASE, "data/raw/WHRadjBMI_Pulit2019_EUR.txt.gz")
HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n"
OUTS = ["CAD", "HF", "T2D", "niHF", "niHFpEF", "niHFrEF"]

def rsid(s):
    return s.split(":", 1)[0]

def part_a_sig(pthr=5e-8):
    out = os.path.join(INSTR, "WHRadjBMI.sig.tsv"); k = 0
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
                    f"{t[ix['Freq_Tested_Allele']]}\t{b}\t{se}\t{p}\t{t[ix['N']]}\tWHRadjBMI\n"); k += 1
    print(f"[WHRadjBMI] wrote {k} genome-wide-sig instruments -> WHRadjBMI.sig.tsv", flush=True)

def part_b_outcomes():
    cf = os.path.join(INSTR, "WHRadjBMI.clumped.tsv")
    if not os.path.exists(cf):
        sys.exit("run 04_clump.R WHRadjBMI first")
    with io.open(cf, encoding="utf-8") as fh:
        fh.readline()
        snps = {ln.split("\t", 1)[0] for ln in fh}
    print(f"[WHRadjBMI] {len(snps)} clumped instruments", flush=True)
    for o in OUTS:
        outf = os.path.join(HARM, f"WHRadjBMI__{o}.outcome.tsv")
        if os.path.exists(outf):
            print(f"  WHRadjBMI->{o}: exists, skip"); continue
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
        print(f"  WHRadjBMI->{o}: {len(rows)}/{len(snps)} SNPs", flush=True)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sig"
    (part_a_sig if mode == "sig" else part_b_outcomes)()
