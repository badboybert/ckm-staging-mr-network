# -*- coding: utf-8 -*-
"""Audit the marker types in the SBP outcome GWAS, and record the answer as an artifact.

Round-2 review item 12 asked for variant-level QC on the disease->SBP edges. Script 62 supplied it
and returned 0 allele mismatches and 0 non-SNV matches across 3,156 matched variants. Reporting that
as "the QC passed" would have been misleading, because part of it is TRUE BY CONSTRUCTION: if the
outcome file contains only SNPs, then an indel or multiallelic mismatch cannot occur no matter how
the extraction was done.

This script measures that directly from the source file rather than assuming it, and writes the
count to results/ so the claim is checkable from the analysis tree afterwards. Before this existed,
the number lived only in a session note -- provable only by re-reading a 136 MB gzip, which is the
same as not being provable.

Out: results/sbp_markertype_audit.txt
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

import io, os, re, sys, gzip, collections

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = P4_BASE
SRC = os.path.join(BASE, "data", "raw", "SBP_Evangelou2018_EUR.txt.gz")
OUT = os.path.join(BASE, "results", "sbp_markertype_audit.txt")

TYPE_RE = re.compile(r":(SNP|INDEL|DEL|INS|MULTI|MNP)\b")

n = 0
kinds = collections.Counter()
with gzip.open(SRC, "rt", errors="replace") as fh:
    header = fh.readline().rstrip("\n").split()
    for line in fh:
        f = line.split()
        if not f:
            continue
        n += 1
        m = TYPE_RE.search(f[0])
        kinds[m.group(1) if m else "untyped"] += 1

snp = kinds.get("SNP", 0)
lines = [
    "=== Marker-type audit of the SBP outcome GWAS (Evangelou et al. 2018) ===",
    "",
    f"source  : {os.path.basename(SRC)}",
    f"columns : {' '.join(header)}",
    f"records : {n:,}",
    "",
    "marker type (from the MarkerName suffix):",
]
for k, v in kinds.most_common():
    lines.append(f"  {k:8s} {v:>12,}   {100.0 * v / n:8.4f}%")
lines += [
    "",
    "WHY THIS IS REPORTED:",
    "  Script 62 found 0 allele mismatches and 0 non-SNV matches among 3,156 matched variants across",
    "  the X->SBP edges. That is a real result, but it is PARTLY STRUCTURAL: with only SNP-typed",
    f"  markers in the outcome file ({snp:,} of {n:,}), an indel or multiallelic mismatch cannot occur",
    "  by construction. With 0 mismatches and 0 strand flips, the 'strict' instrument set therefore",
    "  differs from the primary one ONLY by dropping the 454 palindromic variants -- it tests",
    "  palindrome exclusion, not allele matching, and the manuscript says so rather than describing",
    "  the QC as passed.",
]
open(OUT, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nwrote {os.path.relpath(OUT, BASE)}")
assert n > 0, "no records read"
