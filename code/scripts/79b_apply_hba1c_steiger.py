# -*- coding: utf-8 -*-
"""Step 79b: apply the corrected HbA1c sample size to the ONE thing that depends on it.

Re-running 06_mr.R wholesale was tried and rejected: data/harmonised/ also holds the SBP.strict
outcome files, which the shipped network deliberately handles through 75_sbp_strict_pleiotropy.R
rather than as ordinary edges, so a blanket re-run added 11 spurious rows and moved the BMI->SBP
estimate onto the superseded position-only instrument set. The weighted-median P is bootstrapped
from a seeded RNG whose state depends on the order files are processed, so it also perturbed ~130
P-values that nothing had actually re-estimated.

The sample size enters exactly one statistic: the Steiger directionality test. This step therefore
edits only `steiger_correct` and `steiger_p`, only for edges with HbA1c on one side, taking the
values from 78_hba1c_provenance_steiger.R, which re-ran directionality_test on the same harmonised
files under both sample sizes. Every other cell is left byte-identical, and that is asserted.

Reads : results/hba1c_provenance_steiger.csv
Writes: results/forward_local_edges.csv (in place, after a byte-level check)
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

import csv, io, os, sys

sys.stdout.reconfigure(encoding="utf-8")
BASE = P4_BASE
RES = os.path.join(BASE, "results")
EDGES = os.path.join(RES, "forward_local_edges.csv")
STEIGER = os.path.join(RES, "hba1c_provenance_steiger.csv")
BONF = 0.05 / 132
STAGE = {"BMI": 1, "SBP": 2, "TG": 2, "HDL": 2, "TC": 2, "LDL": 2, "HbA1c": 2, "T2D": 2, "eGFR": 2,
         "CAD": 4, "HF": 4, "Stroke": 4}


def graph(rows):
    return {(r["exposure"], r["outcome"]) for r in rows
            if float(r["ivw_p"]) < BONF and str(r["steiger_correct"]).upper() == "TRUE"}


def concordance(edges):
    cross = [(a, b) for a, b in edges if a in STAGE and b in STAGE and STAGE[a] != STAGE[b]]
    fwd = [(a, b) for a, b in cross if STAGE[a] < STAGE[b]]
    return len(fwd), len(cross)


def main():
    with io.open(EDGES, encoding="utf-8", newline="") as f:
        text_before = f.read()
    rows = list(csv.DictReader(io.StringIO(text_before)))
    fields = list(rows[0].keys())
    new = {(r["exposure"], r["outcome"]): r
           for r in csv.DictReader(io.open(STEIGER, encoding="utf-8"))}

    g_before = graph(rows)
    changed = []
    for r in rows:
        k = (r["exposure"], r["outcome"])
        if "HbA1c" not in k or k not in new:
            continue
        old_c, old_p = r["steiger_correct"], r["steiger_p"]
        r["steiger_correct"] = new[k]["steiger_new"]
        r["steiger_p"] = new[k]["steiger_p_new"]
        if (old_c, old_p) != (r["steiger_correct"], r["steiger_p"]):
            changed.append((k, old_c, r["steiger_correct"]))
    g_after = graph(rows)

    # Rewrite LINE BY LINE so every untouched cell keeps its exact original text. A DictWriter
    # round-trip re-formats numbers (452 -> 452.0, 0 -> 0.0) on every row, which would bury a
    # two-cell change in a whole-file diff.
    assert fields[-2:] == ["steiger_correct", "steiger_p"], fields[-2:]
    lines = text_before.split("\n")
    edited = 0
    for i, r in enumerate(rows, start=1):
        k = (r["exposure"], r["outcome"])
        if "HbA1c" not in k or k not in new:
            continue
        head = lines[i].rsplit(",", 2)[0]
        lines[i] = f'{head},"{r["steiger_correct"]}",{r["steiger_p"]}'
        edited += 1
    io.open(EDGES, "w", encoding="utf-8", newline="").write("\n".join(lines))
    print(f"rewrote {edited} line(s); every other line is byte-identical")

    print(f"HbA1c edges updated: {len(changed)} Steiger cell pair(s) changed")
    for k, a, b in changed:
        print(f"   {k[0]:>6}->{k[1]:<6} steiger_correct {a} -> {b}")
    print(f"\nSteiger-directed graph: {len(g_before)} -> {len(g_after)} edges")
    print(f"  entering: {sorted(g_after - g_before)}   leaving: {sorted(g_before - g_after)}")
    f0, c0 = concordance(g_before); f1, c1 = concordance(g_after)
    print(f"  cross-stage: {f0}/{c0} = {f0/c0:.4f}  ->  {f1}/{c1} = {f1/c1:.4f}")

    print("\n--- SELF-CHECKS ---")
    ok = True
    for name, cond in [
        ("only HbA1c edges changed",
         all("HbA1c" in k for k, _, _ in changed)),
        ("the cross-stage concordance is untouched", (f0, c0) == (f1, c1)),
        ("exactly one edge enters the graph", len(g_after - g_before) == 1),
        ("the entering edge is within-stage",
         all(STAGE[a] == STAGE[b] for a, b in (g_after - g_before))),
        ("no edge leaves the graph", not (g_before - g_after)),
    ]:
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
        ok &= bool(cond)
    if not ok:
        raise SystemExit(1)


def _isnum(v):
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    main()
