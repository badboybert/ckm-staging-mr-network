# -*- coding: utf-8 -*-
"""Step 35a (UKB-free OUTCOME sensitivity): extract UKB-free outcome-GWAS rows matching
each EUR exposure's clumped instrument rsIDs. Streams each big outcome file ONCE.
Providers (all in data/uploads/eur_noukb/):
  CAD    : CARDIoGRAM (Nikpay 2015 VCF, per-SNP SS) + FinnGen R12 CHD (pheweb)
  HF     : HERMES noUKB (METAL, per-SNP TotalSampleSize) + FinnGen R12 (pheweb)
  T2D    : Mahajan noUKB (TSV, n_const) + FinnGen R12 (pheweb)
  Stroke : MEGASTROKE (VCF, n_const)
  eGFR   : CKDGen 2016 (VCF, n_const)
Writes data/harmonised_noukb/<EXP>__<O>_<PROVIDER>.outcome.tsv (schema == step 05).
Checkpoint: a provider is skipped if ALL its per-exposure output files already exist.
N.B. exposure side (Yengo-BMI, Evangelou-SBP) still contains UKB -> this is OUTCOME-side only.
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

import os, io, gzip, math

BASE  = P4_BASE
INSTR = os.path.join(BASE, "data/instruments")
HARM  = os.path.join(BASE, "data/harmonised_noukb")
UP    = os.path.join(BASE, "data/uploads/eur_noukb")
os.makedirs(HARM, exist_ok=True)

EXPS_ALL = ["BMI","SBP","TG","HDL","TC","LDL","HbA1c","CAD","HF","T2D","eGFR"]

def exps_for(outcome):
    return [e for e in EXPS_ALL if e != outcome]

# N constants (affect Steiger directionality ONLY, not beta/P) -- all flagged in output.
N_FINNGEN = 500348   # FinnGen R12 data-freeze total (cohort total; endpoint-specific effective N unknown)
N_MAHAJAN = 231436   # Mahajan 2018 Nat Genet EUR (74,124 cases + 157,312 controls), no-UKB DIAGRAM
N_MEGA    = 446696   # MEGASTROKE all-stroke EUR (non-UKB release)
N_CKDGEN  = 567460   # CKDGen 2016 (Pattaro) eGFRcrea EUR

PROVIDERS = [
  dict(outcome="CAD",    tag="CARDIoGRAM", path=os.path.join(UP,"CAD_CARDIoGRAM_ieu-a-7.vcf.gz"),          type="vcf",     n_const=None),
  dict(outcome="CAD",    tag="FinnGen",    path=os.path.join(UP,"CAD_FinnGenR12_CHD_hg19.gz"),             type="pheweb",  n_const=N_FINNGEN),
  dict(outcome="HF",     tag="HERMES",     path=os.path.join(UP,"HF_HERMES_noUKB.tsv.gz"),                 type="metal",   n_const=None),
  dict(outcome="HF",     tag="FinnGen",    path=os.path.join(UP,"HF_FinnGenR12_hg19.gz"),                  type="pheweb",  n_const=N_FINNGEN),
  dict(outcome="T2D",    tag="Mahajan",    path=os.path.join(UP,"T2D_Mahajan_noUKB_rsid.txt"),             type="mahajan", n_const=N_MAHAJAN),
  dict(outcome="T2D",    tag="FinnGen",    path=os.path.join(UP,"T2D_FinnGenR12_hg19.gz"),                 type="pheweb",  n_const=N_FINNGEN),
  dict(outcome="Stroke", tag="MEGASTROKE", path=os.path.join(UP,"Stroke_MEGASTROKE_ebi-a-GCST005838.vcf.gz"), type="vcf", n_const=N_MEGA),
  dict(outcome="eGFR",   tag="CKDGen",     path=os.path.join(UP,"eGFR_CKDGen_ebi-a-GCST003372.vcf.gz"),    type="vcf",     n_const=N_CKDGEN),
]

def _open(p):
    return gzip.open(p, "rt", errors="replace") if p.endswith(".gz") else open(p, "rt", errors="replace")

def load_clumped_snps(node):
    f = os.path.join(INSTR, f"{node}.clumped.tsv")
    if not os.path.exists(f): return set()
    with io.open(f, encoding="utf-8") as fh:
        fh.readline()
        return {ln.split("\t")[0] for ln in fh}

def _p_from_lp(lp):
    try:
        v = 10.0 ** (-float(lp))
    except (ValueError, OverflowError):
        return "NA"
    if v <= 0.0: v = 1e-300
    return repr(v)

def _std(ea, oa, eaf, beta, se, p, N):
    """validate + standardize one record -> tuple or None"""
    try:
        b = float(beta); s = float(se)
    except (ValueError, TypeError):
        return None
    if s <= 0 or math.isnan(b) or math.isnan(s):
        return None
    ea = ea.strip().upper(); oa = oa.strip().upper()
    if ea in ("", ".") or oa in ("", "."):
        return None
    return (ea, oa, str(eaf), repr(b), repr(s), str(p), str(N))

def collect(provider, need):
    """Stream provider file once; return {snp: std-tuple} for rsIDs in `need`."""
    rows = {}
    t = provider["type"]; n_const = provider["n_const"]
    with _open(provider["path"]) as f:
        if t == "vcf":
            # skip ## meta; header line starts with #CHROM
            for line in f:
                if line.startswith("##"): continue
                if line.startswith("#CHROM"): break
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 10: continue
                snp = c[2]
                if snp not in need: continue
                oa = c[3]; ea = c[4]                 # REF=other, ALT=effect (ES is for ALT)
                fmt = c[8].split(":"); smp = c[9].split(":")
                d = dict(zip(fmt, smp))
                beta = d.get("ES"); se = d.get("SE"); eaf = d.get("AF", "NA")
                p = _p_from_lp(d.get("LP", "nan"))
                N = d.get("SS") if ("SS" in d and d.get("SS") not in (None, "", ".")) else (str(n_const) if n_const else "NA")
                r = _std(ea, oa, eaf, beta, se, p, N)
                if r: rows[snp] = r
        elif t == "pheweb":
            hdr = f.readline().rstrip("\n").split("\t")   # #chrom pos ref alt rsids ... pval mlogp beta sebeta af_alt
            ix = {k: hdr.index(k) for k in ("ref","alt","rsids","pval","beta","sebeta","af_alt")}
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) <= ix["af_alt"]: continue
                rs = c[ix["rsids"]].split(",")[0]
                if not rs.startswith("rs") or rs not in need: continue
                r = _std(c[ix["alt"]], c[ix["ref"]], c[ix["af_alt"]], c[ix["beta"]], c[ix["sebeta"]], c[ix["pval"]], n_const)
                if r: rows[rs] = r
        elif t == "metal":
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {k: hdr.index(k) for k in ("Allele1","Allele2","Freq1","Effect","StdErr","P.value","TotalSampleSize","final_rsid")}
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) <= ix["final_rsid"]: continue
                rs = c[ix["final_rsid"]]
                if rs not in need: continue
                r = _std(c[ix["Allele1"]], c[ix["Allele2"]], c[ix["Freq1"]], c[ix["Effect"]], c[ix["StdErr"]], c[ix["P.value"]], c[ix["TotalSampleSize"]])
                if r: rows[rs] = r
        elif t == "mahajan":
            hdr = f.readline().rstrip("\n").split("\t")    # rsID BP SNP CHR EA NEA EAF Beta SE Pvalue
            ix = {k: hdr.index(k) for k in ("rsID","EA","NEA","EAF","Beta","SE","Pvalue")}
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) <= ix["Pvalue"]: continue
                rs = c[ix["rsID"]]
                if not rs.startswith("rs") or rs not in need: continue
                r = _std(c[ix["EA"]], c[ix["NEA"]], c[ix["EAF"]], c[ix["Beta"]], c[ix["SE"]], c[ix["Pvalue"]], n_const)
                if r: rows[rs] = r
    return rows

def outfile(exp, outcome, tag):
    return os.path.join(HARM, f"{exp}__{outcome}_{tag}.outcome.tsv")

for prov in PROVIDERS:
    o = prov["outcome"]; tag = prov["tag"]
    exps = exps_for(o)
    expected = [outfile(e, o, tag) for e in exps]
    if all(os.path.exists(p) and os.path.getsize(p) > 0 for p in expected):
        print(f"[skip] {o}/{tag}: all {len(expected)} outcome files present", flush=True)
        continue
    if not os.path.exists(prov["path"]):
        print(f"[MISS] {o}/{tag}: file not found {prov['path']} -- skipping", flush=True)
        continue
    # union of needed SNPs across all exposures targeting this outcome
    need = set()
    exp_snps = {}
    for e in exps:
        s = load_clumped_snps(e); exp_snps[e] = s; need |= s
    print(f"[stream] {o}/{tag}: need {len(need)} unique instrument SNPs across {len(exps)} exposures ...", flush=True)
    rows = collect(prov, need)
    print(f"[stream] {o}/{tag}: matched {len(rows)}/{len(need)} SNPs", flush=True)
    for e in exps:
        outf = outfile(e, o, tag); m = 0
        with io.open(outf, "w", encoding="utf-8") as w:
            w.write("SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tpval\tN\tPhenotype\n")
            for s in exp_snps[e]:
                if s in rows:
                    ea, oa, eaf, b, se, p, N = rows[s]
                    w.write(f"{s}\t{ea}\t{oa}\t{eaf}\t{b}\t{se}\t{p}\t{N}\t{o}\n"); m += 1
        print(f"   {e}->{o} [{tag}]: {m} outcome SNPs", flush=True)
print("done")
