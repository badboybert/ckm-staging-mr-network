# -*- coding: utf-8 -*-
"""Step 1: exhaustive local GWAS inventory + CKM-node coverage + gap list.
Scans the whole project for GWAS-like summary-stat files, classifies by trait node
and ancestry, parses the header + first data row to identify columns. Observable
(prints progress) and writes durable CSVs to manifest/."""
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

import os, gzip, io, csv, re

ROOT = CKM_ROOT
OUT  = os.path.join(P4_BASE, "manifest")

# --- CKM node trait keyword patterns. Matched against TOKENS (basename split on
# non-alphanumeric), so short codes like 'hf','tg','cad' match 'HF_GCST...','tg.gz'. ---
# Each entry: (node, set of exact-token codes, list of substring patterns).
NODE_PATTS = [
    ("BMI",    {"bmi"}, ["body_mass", "bodymass"]),
    ("SBP",    {"sbp","dbp","bp","icbp"}, ["systolic", "blood_pressure", "bloodpressure", "hypertension", "htn"]),
    ("TG",     {"tg","trig"}, ["triglyc"]),
    ("HDL",    {"hdl"}, []),
    ("LDL",    {"ldl"}, []),
    ("TC",     {"tc","chol"}, ["total_chol", "totalchol", "cholesterol"]),
    ("HbA1c",  {"hba1c","a1c"}, ["glycated"]),
    ("FG_FI",  {"fg","fi"}, ["fasting_gluc","fasting_insul","fastinggluc","fastinginsul","glucose","insulin","homa"]),
    ("T2D",    {"t2d","t2dm","diagram"}, ["type2","type_2","diabet","mahajan","agen_t2d","spracklen"]),
    ("eGFR",   {"egfr","gfr"}, ["creatinine","ckdgen"]),
    ("CKD",    {"ckd"}, ["chronic_kidney","kidney_disease"]),
    ("UACR",   {"uacr"}, ["albumin","albuminuria"]),
    ("NAFLD",  {"nafld","masld"}, ["liver_fat","liverfat","steato","pnpla3"]),
    ("ALT",    {"alt"}, ["alanine"]),
    ("AST",    {"ast"}, ["aspartate"]),
    ("GGT",    {"ggt"}, ["gamma_glut","gammaglut"]),
    ("CAD",    {"cad","chd","mi"}, ["coronary","cardiogram","myocard","ischaemic_heart","ischemic_heart"]),
    ("HF",     {"hf"}, ["heart_fail","heartfail","hermes"]),
    ("Stroke", {"stroke","is"}, ["megastroke","ischaemic_stroke","ischemic_stroke"]),
    ("AF",     {"af","afib"}, ["atrial_fib","atrialfib"]),
]
ANC_TOKENS_EAS = {"eas","bbj","japan","japanese","agen","korean","tpmi","taiwan","kadoorie","tommo","tomo","asian","han"}
ANC_TOKENS_EUR = {"eur","european","ukb","finn","finngen","giant","glgc","icbp","cardiogram","hermes","ckdgen","megastroke","hunt"}
# QTL / non-GWAS-node markers — matched against BASENAME TOKENS only (path-independent)
EXCLUDE_TOK = {"eqtl","mqtl","pqtl","meqtl","methyl","betas","soft","godmc","susztak","aptamer","seqid","coloc","betasgeo"}
EXCLUDE_RE  = r"(ensg\d|^oid\d|qtd\d|_family|series_matrix|__enst|leads)"

SUMSTAT_EXT = (".tsv.gz", ".txt.gz", ".vcf.gz", ".gz", ".sumstats.gz", ".tsv", ".txt")

def tokenize(s):
    return set(t for t in re.split(r"[^a-z0-9]+", s.lower()) if t)

def classify(path):
    low = path.replace("\\", "/").lower()
    base = os.path.basename(low)
    btok = tokenize(base)
    # exclude molecular-QTL/omics files by basename tokens
    if btok & EXCLUDE_TOK or re.search(EXCLUDE_RE, base):
        return None, None
    node = None
    for n, codes, subs in NODE_PATTS:
        if (btok & codes) or any(sub in base for sub in subs):
            node = n; break
    # ancestry: prefer explicit token anywhere in path
    ptok = tokenize(low)
    if ptok & ANC_TOKENS_EAS:
        anc = "EAS"
    elif ptok & ANC_TOKENS_EUR:
        anc = "EUR"
    else:
        anc = "EUR?"
    return node, anc

def peek_header(path, nlines=2):
    try:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt", errors="replace") as f:
            rows = [f.readline().strip() for _ in range(nlines)]
        return rows
    except Exception as e:
        return [f"<read-error: {e}>"]

def guess_cols(header):
    h = header.lower()
    def has(*ks): return any(k in h for k in ks)
    return {
        "chr":  has("chr", "chrom", "#chr"),
        "pos":  has("pos", "bp", "base_pair", "position"),
        "rsid": has("rsid", "snp", "variant", "rs_"),
        "ea":   has("effect_allele", "alt", "a1", "effectallele", "ea"),
        "beta": has("beta", "effect", "b\t", "\tb", "or"),
        "se":   has("se", "standard_error", "standarderror"),
        "p":    has("p_value", "pval", "\tp\t", "p-value", "min_log10", "logp"),
        "n":    has("\tn", "n_total", "samplesize", "\tn_", "neff"),
    }

rows = []
scanned = 0
print("Scanning", ROOT, "...", flush=True)
for dirpath, _, files in os.walk(ROOT):
    for fn in files:
        if not fn.lower().endswith(SUMSTAT_EXT):
            continue
        full = os.path.join(dirpath, fn)
        scanned += 1
        node, anc = classify(full)          # returns (None, None) for excluded/unmatched
        if node is None:
            continue
        try:
            size = os.path.getsize(full)
        except OSError:
            size = -1
        if size < 20_000_000:   # <20MB: cis-region extract / coloc clutter, not a full GWAS
            continue
        hdr = peek_header(full)
        cols = guess_cols(hdr[0]) if hdr else {}
        ncol_ok = sum(bool(v) for v in cols.values())
        rows.append(dict(
            node=node, ancestry=anc, size_mb=round(size/1e6, 1),
            path=full.replace("\\", "/"),
            rel=os.path.relpath(full, ROOT).replace("\\", "/"),
            header=hdr[0][:200] if hdr else "",
            col_score=ncol_ok,
            has_beta=cols.get("beta"), has_se=cols.get("se"),
            has_p=cols.get("p"), has_rsid=cols.get("rsid"),
        ))
        print(f"  [{node:6s}|{anc:4s}] {round(size/1e6):5d}MB  {os.path.relpath(full, ROOT)[:80]}", flush=True)

print(f"\nScanned {scanned} candidate files; {len(rows)} classified as CKM-node GWAS.", flush=True)

# write full manifest
os.makedirs(OUT, exist_ok=True)
with io.open(os.path.join(OUT, "local_gwas_manifest.csv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                       ["node","ancestry","size_mb","path","rel","header","col_score","has_beta","has_se","has_p","has_rsid"])
    w.writeheader()
    for r in sorted(rows, key=lambda r: (r["node"], r["ancestry"], -r["size_mb"])):
        w.writerow(r)

# coverage matrix: node x ancestry -> best (largest, col-complete) file
NODES = [n for n, *_ in NODE_PATTS]
cov = {}
for r in rows:
    key = (r["node"], "EAS" if r["ancestry"] == "EAS" else "EUR")
    cur = cov.get(key)
    if cur is None or (r["col_score"], r["size_mb"]) > (cur["col_score"], cur["size_mb"]):
        cov[key] = r

with io.open(os.path.join(OUT, "node_coverage.csv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["node", "EUR_local", "EUR_file", "EAS_local", "EAS_file"])
    for n in NODES:
        e = cov.get((n, "EUR")); a = cov.get((n, "EAS"))
        w.writerow([n,
                    "YES" if e else "no", e["rel"] if e else "",
                    "YES" if a else "no", a["rel"] if a else ""])

print("\n=== CKM NODE COVERAGE (local) ===")
print(f"{'node':7s} {'EUR':4s} {'EAS':4s}  best-EUR-file")
for n in NODES:
    e = cov.get((n, "EUR")); a = cov.get((n, "EAS"))
    print(f"{n:7s} {'YES' if e else '  -':4s} {'YES' if a else '  -':4s}  {(e['rel'][:70] if e else '')}")
print("\nManifest ->", os.path.join(OUT, "local_gwas_manifest.csv"))
print("Coverage ->", os.path.join(OUT, "node_coverage.csv"))
