# -*- coding: utf-8 -*-
"""Step 59: genome-wide-strict staging robustness variant (replaces the stale staging_forward.txt).

The earlier "genome-wide-strict FORWARD" variant restricted the graph to forward-direction edges and
then reported concordance 1.000 — concordance was 1.000 by construction, a circular test, and the file
(results/staging_forward.txt, 2026-07-11) was built from the pre-completion 51-edge forward network.

This replaces it with an honest, non-circular robustness check on the CURRENT bidirectional network:
raise the edge-inclusion threshold from network Bonferroni (0.05/132) to genome-wide significance
(P < 5e-8), keep Steiger-directed staged edges, and score cross-stage concordance with an EXACT
permutation P (enumerated, matching 07_staging.py). At the stricter bar one backward edge (the
reverse CAD->SBP index-event edge, P = 7e-10) survives while CAD->T2D (P = 1.6e-6) drops, so the
concordance is 1.000 (19/19) since the strict-SBP promotion: the one backward edge, CAD->SBP,
falls just above the genome-wide threshold under strict allele matching (P = 1.11e-07 vs 5e-8) and
drops out. The rise from the primary 0.926 is a threshold effect, not independent evidence.

Writes results/staging_forward.txt (the name is kept so fig1.R's Figure-1d panel reads it unchanged;
the panel label is updated from "strict (forward)" to "strict"). Deterministic; no randomness.
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

import csv, io, os
from collections import Counter

BASE = P4_BASE
GW = 5e-8
STAGE = {"BMI": 1, "SBP": 2, "TG": 2, "HDL": 2, "TC": 2, "LDL": 2, "HbA1c": 2, "T2D": 2, "eGFR": 2,
         "CAD": 4, "HF": 4, "Stroke": 4}


def load(p):
    return list(csv.DictReader(io.open(p, encoding="utf-8")))


def sig_gw(rows):
    out = []
    for r in rows:
        try:
            p = float(r["ivw_p"])
        except (ValueError, KeyError):
            continue
        if p < GW and str(r["steiger_correct"]).upper() == "TRUE" \
           and r["exposure"] in STAGE and r["outcome"] in STAGE:
            out.append((r["exposure"], r["outcome"]))
    return out


def concordance(edges):
    cross = [(a, b) for a, b in edges if STAGE[a] != STAGE[b]]
    fwd = [(a, b) for a, b in cross if STAGE[a] < STAGE[b]]
    back = [(a, b) for a, b in cross if STAGE[a] > STAGE[b]]
    return (len(fwd) / len(cross) if cross else float("nan")), fwd, back, cross


def _msperm(counts):
    n = sum(counts.values()); cur = [None] * n
    def rec(i):
        if i == n:
            yield tuple(cur); return
        for v in sorted(counts):
            if counts[v]:
                counts[v] -= 1; cur[i] = v
                yield from rec(i + 1)
                counts[v] += 1
    return rec(0)


def exact_perm_p(edges):
    nodes = sorted({x for e in edges for x in e})
    ob, _, _, _ = concordance(edges)
    ge = tot = 0
    for pm in _msperm(dict(Counter(STAGE[n] for n in nodes))):
        tot += 1
        smap = dict(zip(nodes, pm))
        cc = [(a, b) for a, b in edges if smap[a] != smap[b]]
        if cc and sum(1 for a, b in cc if smap[a] < smap[b]) / len(cc) >= ob - 1e-12:
            ge += 1
    return ge / tot, ge, tot, len(nodes)


def main():
    rows = load(os.path.join(BASE, "results/forward_local_edges.csv"))
    sig = sig_gw(rows)
    conc, fwd, back, cross = concordance(sig)
    p, ge, tot, nn = exact_perm_p(sig)

    L = []
    def out(s=""):
        L.append(s); print(s, flush=True)

    out("=== Genome-wide-strict staging robustness variant (on the completed 132-edge network) ===")
    out("Edge-inclusion threshold raised from network Bonferroni (0.05/132) to genome-wide "
        "significance (P < 5e-8); Steiger-directed staged edges only.")
    out(f"Significant Steiger-correct causal edges: {len(sig)} over {nn} nodes")
    out("")
    out(f"[3] Cross-stage edges: {len(cross)}  | forward(low->high AHA)={len(fwd)}  "
        f"backward(high->low)={len(back)}")
    out(f"    Stage-order CONCORDANCE = {conc:.3f}")
    out(f"    Backward (AHA-discordant) edges: {back if back else 'none'}")
    out("    (Since X->SBP was promoted to strict allele-compatible matching, the reverse CAD->SBP "
        "index-event edge no longer clears the genome-wide bar: P = 1.11e-07 against 5e-8, so it "
        "drops out along with CAD->T2D at P = 1.6e-6. The resulting 1.000 is HIGHER than the primary "
        "0.926 only because the discordant edge fell below a threshold, and is not independent "
        "evidence for the ordering.)")
    out("")
    # Full-precision P first, so fig1.R's `empirical p = ([0-9.eE+-]+)` regex captures the value (not
    # the numerator of a leading fraction) and its %.1e render (7.6e-04) matches the prose.
    out(f"[5] Concordance vs label-permutation null (EXACT enumeration of {tot:,} distinct labelings): "
        f"empirical p = {p:.6f}  (exact, {ge}/{tot})")
    out(f"    => {'concordant with AHA order beyond chance' if p < 0.05 else 'NOT beyond chance'}")
    out("")
    out("NOTE: tests necessary POPULATION-LEVEL causal-ordering conditions, not within-person progression.")

    io.open(os.path.join(BASE, "results/staging_forward.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")

    # self-checks against the adopted (Option 1) targets
    print("\n--- SELF-CHECKS ---", flush=True)
    ok = True
    for name, cond, got in [
        ("concordance 1.000 (19/19) after the strict-SBP promotion",
         f"{conc:.3f}" == "1.000" and (len(fwd), len(cross)) == (19, 19),
         (f"{conc:.3f}", len(fwd), len(cross))),
        ("no backward edge clears the genome-wide bar under strict matching", back == [], back),
        ("exact P = 1/1,320 = 7.6e-4", (ge, tot) == (1, 1320) and f"{p:.1e}" == "7.6e-04",
         (ge, tot, f"{p:.1e}")),
    ]:
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}  (got {got})", flush=True)
        ok &= bool(cond)
    print(f"\nALL TARGETS REPRODUCED: {ok}", flush=True)
    print("wrote results/staging_forward.txt", flush=True)


if __name__ == "__main__":
    main()
