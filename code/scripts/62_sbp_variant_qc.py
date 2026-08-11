# -*- coding: utf-8 -*-
"""Step 62: variant-level QC for the disease->SBP edges, and a strict-matched re-extraction.

Why this is load-bearing. CAD->SBP is the only backward cross-stage edge that survives the
genome-wide-strict threshold (concordance 0.950), one of only two backward edges in the primary
network, and the single backward edge in the UK-Biobank-free MAIN network. Every X->SBP outcome was
extracted by POSITION ONLY:

    37_complete_network.py: rsID -> (chr,pos) via the build37 bim, then SBP rows matched on
    (chr,pos) with `wanted.setdefault(cp, rs)` and `if rs in sbp_at: continue`.

Three silent failure modes follow, none of which the extractor checks:
  (a) NO ALLELE CHECK. The SBP record's Allele1/Allele2 at that position are taken as the outcome
      alleles whatever they are; a different variant at the same coordinate is accepted.
  (b) The Evangelou MarkerName is `chr:pos:type`; the TYPE field is discarded, so an indel or a
      multiallelic record can be matched to an SNV instrument.
  (c) Position collisions (two instrument rsIDs at one coordinate) and multiallelic SBP records
      (several rows at one coordinate) are both resolved by "keep the first seen".

`harmonise_data(action=2)` downstream does align alleles and drop ambiguous palindromes, so this is
a bound on the risk rather than proof of error. This script quantifies it: it classifies every
matched variant, then writes a STRICT outcome file keeping only variants whose alleles are
unambiguously concordant with the exposure, for re-estimation by 63_sbp_strict_mr.R.

Out: results/sbp_variant_qc.{txt,csv} and data/harmonised/*__SBP.strict.outcome.tsv
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

import csv, gzip, io, os, collections

BASE = P4_BASE
ROOT = CKM_ROOT
INSTR = os.path.join(BASE, "data/instruments")
HARM = os.path.join(BASE, "data/harmonised")
RES = os.path.join(BASE, "results")
BIM = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur.bim")
SBP_GWAS = os.path.join(BASE, "data/raw/SBP_Evangelou2018_EUR.txt.gz")

SBP_EXPS = ["BMI", "TG", "HDL", "TC", "LDL", "HbA1c", "CAD", "HF", "T2D", "eGFR", "Stroke"]
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def read_tsv(p):
    with io.open(p, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def is_snv(a):
    return len(a) == 1 and a in COMP


def classify(e_ea, e_oa, o_ea, o_oa):
    """Allele relationship between the exposure instrument and the SBP record at that position."""
    if not (is_snv(e_ea) and is_snv(e_oa)):
        return "exposure_not_snv"
    if not (is_snv(o_ea) and is_snv(o_oa)):
        return "outcome_not_snv"          # indel / multi-character allele at the same coordinate
    e, o = {e_ea, e_oa}, {o_ea, o_oa}
    if e == o:
        return "palindromic" if (e_ea, e_oa) in PALINDROMIC else "match"
    if {COMP[e_ea], COMP[e_oa]} == o:
        return "palindromic" if (e_ea, e_oa) in PALINDROMIC else "strand_flip"
    return "allele_mismatch"              # a DIFFERENT variant at the same coordinate


def main():
    # ---- instrument rsID -> position, and detect position collisions ------------------------
    union = set()
    exp_snps = {}
    for e in SBP_EXPS:
        f = os.path.join(INSTR, f"{e}.clumped.tsv")
        rows = read_tsv(f)
        exp_snps[e] = {r["SNP"]: r for r in rows}
        union |= set(exp_snps[e])
    pos_of, at_pos = {}, collections.defaultdict(list)
    with io.open(BIM, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.split()
            if p[1] in union:
                pos_of.setdefault(p[1], (p[0], p[3]))
                at_pos[(p[0], p[3])].append(p[1])
    collisions = {k: v for k, v in at_pos.items() if len(v) > 1}

    # ---- multiallelic / typed SBP records at the wanted positions ----------------------------
    wanted = {pos_of[rs] for rs in pos_of}
    sbp_rows = collections.defaultdict(list)
    with gzip.open(SBP_GWAS, "rt", errors="replace") as f:
        hdr = f.readline().split()
        ix = {c: hdr.index(c) for c in ["MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P"]}
        for line in f:
            t = line.split()
            try:
                mk = t[ix["MarkerName"]].split(":")
            except IndexError:
                continue
            if len(mk) < 2:
                continue
            key = (mk[0], mk[1])
            if key in wanted:
                sbp_rows[key].append((t[ix["MarkerName"]], t[ix["Allele1"]].upper(),
                                      t[ix["Allele2"]].upper(), t[ix["Freq1"]],
                                      t[ix["Effect"]], t[ix["StdErr"]], t[ix["P"]]))
    multi = {k: v for k, v in sbp_rows.items() if len(v) > 1}

    # ---- classify every shipped X->SBP outcome variant ---------------------------------------
    L, per_edge = [], {}
    def out(s=""):
        L.append(s); print(s, flush=True)

    out("=== Variant-level QC for the disease->SBP edges ===")
    out("")
    out("EXTRACTION AUDIT (37_complete_network.py, X->SBP block):")
    out(f"  instrument rsIDs mapped to build37 positions : {len(pos_of)}")
    out(f"  positions carrying >1 instrument rsID (collisions, first-seen wins) : {len(collisions)}")
    if collisions:
        for k, v in sorted(collisions.items())[:10]:
            out(f"      chr{k[0]}:{k[1]} -> {v}")
    out(f"  wanted positions with >1 SBP record (multiallelic, first-seen wins) : {len(multi)}")
    if multi:
        for k, v in sorted(multi.items())[:10]:
            out(f"      chr{k[0]}:{k[1]} -> {[m[0] for m in v]}")
    out("")
    out("PER-EDGE ALLELE CONCORDANCE (exposure instrument vs the SBP record taken at its position):")
    out(f"{'edge':14s} {'n':>5s} {'match':>6s} {'flip':>5s} {'palin':>6s} {'mismatch':>9s} {'notSNV':>7s}")
    rows_csv = []
    for e in SBP_EXPS:
        p = os.path.join(HARM, f"{e}__SBP.outcome.tsv")
        if not os.path.exists(p):
            continue
        cnt = collections.Counter()
        keep = []
        for r in read_tsv(p):
            src = exp_snps[e].get(r["SNP"])
            if src is None:
                cnt["no_exposure_row"] += 1
                continue
            c = classify(src["effect_allele"].upper(), src["other_allele"].upper(),
                         r["effect_allele"].upper(), r["other_allele"].upper())
            cnt[c] += 1
            if c in ("match", "strand_flip"):
                keep.append(r)
            rows_csv.append(dict(edge=f"{e}->SBP", SNP=r["SNP"],
                                 exp_ea=src["effect_allele"], exp_oa=src["other_allele"],
                                 out_ea=r["effect_allele"], out_oa=r["other_allele"],
                                 classification=c,
                                 pos_collision=pos_of.get(r["SNP"]) in collisions,
                                 multiallelic=pos_of.get(r["SNP"]) in multi))
        n = sum(cnt.values())
        per_edge[e] = (cnt, n, len(keep))
        out(f"{e+'->SBP':14s} {n:5d} {cnt['match']:6d} {cnt['strand_flip']:5d} "
            f"{cnt['palindromic']:6d} {cnt['allele_mismatch']:9d} "
            f"{cnt['outcome_not_snv']+cnt['exposure_not_snv']:7d}")
        # strict file: unambiguous concordant variants only
        sp = os.path.join(HARM, f"{e}__SBP.strict.outcome.tsv")
        with io.open(sp, "w", encoding="utf-8") as w:
            w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
            for r in keep:
                w.write("\t".join([r["SNP"], r["effect_allele"], r["other_allele"], r["eaf"],
                                   r["beta"], r["se"], r["pval"], r["N"], "SBP"]) + "\n")

    tot = collections.Counter()
    for cnt, _, _ in per_edge.values():
        tot.update(cnt)
    out("")
    out(f"TOTAL across all X->SBP: {sum(tot.values())} matched variants — "
        f"match {tot['match']}, strand-flip {tot['strand_flip']}, palindromic {tot['palindromic']}, "
        f"allele-mismatch {tot['allele_mismatch']}, non-SNV {tot['outcome_not_snv']+tot['exposure_not_snv']}")
    out("")
    out("READ:")
    out("- 'allele_mismatch' and 'non-SNV' are variants where the SBP record at that coordinate is a")
    out("  DIFFERENT variant from the instrument. Position-only matching cannot detect these; they are")
    out("  the exposure the reviewer identified.")
    out("- 'palindromic' variants are strand-ambiguous. harmonise_data(action=2) drops those with")
    out("  intermediate frequency downstream, so they are a bound rather than an error.")
    out("- The strict files keep only 'match' and 'strand_flip'. 63_sbp_strict_mr.R re-estimates")
    out("  every X->SBP edge on them, CAD->SBP above all.")

    io.open(os.path.join(RES, "sbp_variant_qc.txt"), "w", encoding="utf-8").write("\n".join(L))
    with io.open(os.path.join(RES, "sbp_variant_qc.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_csv[0].keys()))
        w.writeheader()
        w.writerows(rows_csv)
    print("\nwrote results/sbp_variant_qc.{txt,csv} and *__SBP.strict.outcome.tsv")


if __name__ == "__main__":
    main()
