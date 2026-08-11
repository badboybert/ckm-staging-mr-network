# -*- coding: utf-8 -*-
"""Step 34 (ITEM 1 correction check): FULL-DENSITY LDL+HDL+TG -> CAD MVMR, no ApoB coverage restriction.
Tests the verifier's claim that HDL RETAINS a significant direct effect (~-0.166, p1.3e-5) at full
instrument density, i.e. the 124-SNP 'HDL dissolves' was a power artifact of the ApoB complete-case funnel.
All of LDL/HDL/TG/CAD are rsID-keyed on the hg19 panel -> NO liftover needed."""
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

import io, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_formats import FORMATS, ROOT, _open

BASE  = os.path.join(ROOT, "paper 4/independent_build")
INSTR = os.path.join(BASE, "data/instruments")
MVDIR = os.path.join(BASE, "data/mvmr"); os.makedirs(MVDIR, exist_ok=True)
PLINK = os.path.join(ROOT, "lhcmr_eas_refs/tools/plink.exe")
PANEL = os.path.join(ROOT, "paper 6/analysis/g1000_eur/g1000_eur")
BIM   = PANEL + ".bim"
EXPO = ["LDL", "HDL", "TG"]; NODES = ["LDL", "HDL", "TG", "CAD"]

def read_clumped(node):
    d = {}
    with io.open(os.path.join(INSTR, node + ".clumped.tsv"), encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t"); ip = h.index("pval"); isnp = h.index("SNP")
        for line in f:
            t = line.rstrip("\n").split("\t")
            try: d[t[isnp]] = min(d.get(t[isnp], 1.0), float(t[ip]))
            except (ValueError, IndexError): continue
    return d

minp = {}
for node in EXPO:
    for rs, p in read_clumped(node).items(): minp[rs] = min(minp.get(rs, 1.0), p)
print(f"[union] {len(minp)} LDL+HDL+TG instruments")
assoc = os.path.join(MVDIR, "lipid_union_assoc.txt")
with io.open(assoc, "w", encoding="utf-8") as w:
    w.write("SNP\tP\n")
    for rs, p in minp.items(): w.write(f"{rs}\t{p}\n")
out_prefix = os.path.join(MVDIR, "lipid_union_clump")
subprocess.run([PLINK, "--bfile", PANEL, "--clump", assoc, "--clump-p1", "1", "--clump-r2", "0.001",
                "--clump-kb", "10000", "--clump-snp-field", "SNP", "--clump-field", "P", "--out", out_prefix],
               capture_output=True, text=True)
clumped = set()
with io.open(out_prefix + ".clumped", encoding="utf-8") as f:
    f.readline()
    for line in f:
        parts = line.split()
        if len(parts) >= 3 and parts[2].startswith("rs"): clumped.add(parts[2])
print(f"[clump] {len(clumped)} LD-independent lipid instruments")

ref = {}
with io.open(BIM, encoding="utf-8", errors="replace") as f:
    for line in f:
        p = line.split()
        if len(p) >= 6 and p[1] in clumped: ref[p[1]] = (p[0], p[3], p[4].upper(), p[5].upper())

def palindromic(a1, a2):
    s = {a1, a2}; return s == {"A", "T"} or s == {"C", "G"}

def extract_standard(node):
    cfg = FORMATS[node]; path = os.path.join(ROOT, cfg["path"]); C = cfg["cols"]
    want = set(ref.keys()); got = {}
    with _open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: hdr.index(v) for k, v in C.items() if v is not None and v in hdr}
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) <= idx.get("snp", 0): t = line.split()
            try: rs = t[idx["snp"]].strip()
            except IndexError: continue
            if rs not in want: continue
            try:
                ea = t[idx["ea"]].upper(); oa = t[idx["oa"]].upper()
                beta = float(t[idx["beta"]]); se = float(t[idx["se"]])
                eaf = float(t[idx["eaf"]]) if "eaf" in idx else float("nan")
            except (ValueError, IndexError): continue
            if se <= 0: continue
            got[rs] = (ea, oa, eaf, beta, se)
            if len(got) == len(want): break
    return got

cache = {n: extract_standard(n) for n in NODES}
for n in NODES: print(f"[extract] {n}: {len(cache[n])}/{len(ref)}")

def aligned(node, rs):
    a1, a2 = ref[rs][2], ref[rs][3]
    rec = cache[node].get(rs)
    if rec is None: return None
    ea, oa, eaf, beta, se = rec
    if {ea, oa} != {a1, a2}: return None
    return (beta, se) if ea == a1 else (-beta, se)

rows = []
for rs, (c, p, a1, a2) in ref.items():
    if palindromic(a1, a2): continue
    vals = {}; ok = True
    for node in NODES:
        a = aligned(node, rs)
        if a is None: ok = False; break
        vals[node] = a
    if ok: rows.append((rs, vals))
out = os.path.join(MVDIR, "model_LDL_HDL_TG_fulldens.tsv")
with io.open(out, "w", encoding="utf-8", newline="") as w:
    w.write("SNP\tbx_LDL\tse_LDL\tbx_HDL\tse_HDL\tbx_TG\tse_TG\tby\tsey\n")
    for rs, v in rows:
        w.write(f"{rs}\t{v['LDL'][0]:.6g}\t{v['LDL'][1]:.6g}\t{v['HDL'][0]:.6g}\t{v['HDL'][1]:.6g}\t"
                f"{v['TG'][0]:.6g}\t{v['TG'][1]:.6g}\t{v['CAD'][0]:.6g}\t{v['CAD'][1]:.6g}\n")
print(f"[model] LDL+HDL+TG->CAD full density -> {len(rows)} SNPs ({os.path.basename(out)})")
